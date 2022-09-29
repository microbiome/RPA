#source("loadData.R")
source("loadData2.R")

require(RPA)
fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}
#cels <- sample(list.celfiles("~/data/RPA/LukkAtlas/", full.names = T), 300)

#emat.rpa <- exprs(rpa(cel.files = cels))
#colnames(emat.rpa) <- sapply(strsplit(colnames(emat.rpa), "/"), function (x) {x[[2]]})
#colnames(emat.rpa) <- sapply(strsplit(colnames(emat.rpa), "\\."), function (x) {x[[1]]})

load(emat.files[["rma"]]) # emat.rma
load("~/Rpackages/RPA/OnlineLearning/emat.rpa.online.storage.RData") 
emat.rpa <- emat.rpa.online.storage

set.seed(342)

N <- 10
classif.method <- "random.forest"
accs <- matrix(NA, nrow = N, ncol = 2)
colnames(accs) <- c("rma", "rpa")
for (cvn in 1:N) {

  # Split training/test data
  train.samples <- sample(colnames(emat.rpa), floor(ncol(emat.rpa)*0.8))
  test.samples  <- setdiff(colnames(emat.rpa), train.samples)

  for (method in c("rpa", "rma")) {

  # Pick data for the selected method
    if (method == "rpa") { emat <- emat.rpa}    # version used in the other experiments
    if (method == "rma") { emat <- emat.rma[, colnames(emat.rpa)]}

  # Select random subset of genes to speed up
  source("~/scripts/R/Affy.R")
	
  all.sets <- removeAFFX(rownames(emat))
  gene.subset <- sample(all.sets, round(0.1*length(all.sets)))
  #gene.subset <- sample(all.sets, 1000)

     # Random forest
     library(randomForest)

     # Training data
     train.data <- data.frame(cbind(t(emat[gene.subset, train.samples])))
     train.data$classInfo <- factor(cl[train.samples])

     # Test data
     test.data <- t(emat[gene.subset, test.samples])
     colnames(test.data) <- paste("X", colnames(test.data), sep = "")

     # Fit the model
     rf <- randomForest(classInfo ~ ., data = train.data, importance = FALSE, proximity = FALSE, ntree = 500)

     # Test prediction on test data
     pred <- predict(rf, test.data)		  
     acc <- mean(pred == cl[test.samples])

     accs[cvn, method] <- acc

  print(paste(cvn, method, acc))

}

}

#save(accs, class.minsize, file = "accs.minsize5.RData")

library(ggplot2)
library(reshape)
o <- order(apply(accs, 2, mean))
accs <- accs[,o]
df <- data.frame(accs); df$fold <- 1:10
df <- melt(df, id = c("fold"))
df$method <- df$variable
ggplot(df) + aes(x = method, y = value, group = method, fill = method) + geom_boxplot()
