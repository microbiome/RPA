source("loadData.R")

set.seed(342)

N <- 10
ntree <- 500
fs <- 0.1

classif.method <- "random.forest"
accs <- matrix(NA, nrow = N, ncol = length(emat.files))
colnames(accs) <- names(emat.files)
for (cvn in 1:N) {

  # Split training/test data
  train.samples <- sample(samples, floor(ncol(emat)*0.8))
  test.samples  <- setdiff(samples, train.samples)

  for (method in names(emat.files)) {
  #for (method in c("rpa.online", "rma", "frma")) {

  # Pick data for the selected method
    load(emat.files[[method]])
    if (method == "rpa") { emat <- emat.rpa[, s]}    # version used in the other experiments
    if (method == "rpa.online") { emat <- emat.rpa.online.storage[, s]} # latest version
    if (method == "rma") { emat <- emat.rma[, s]}
    if (method == "frma") { emat <- emat.frma[, s]}
    if (method == "frmab") { emat <- emat.frma.batch[, s]}
    if (method == "mas") { emat <- emat.MAS5[, s]}

  # Select random subset of genes to speed up
  source("~/scripts/R/Affy.R")
	
  all.sets <- removeAFFX(rownames(emat))
  gene.subset <- sample(all.sets, round(fs*length(all.sets)))
  #gene.subset <- sample(all.sets, 1000)

  if (classif.method == "knn") {

    # KNN / TODO
    library(class)
    res <- knn(t(emat[, train.samples]), t(emat[, test.samples]), factor(cl[train.samples]), k = 1, prob = FALSE)
    attributes(.Last.value)

  } else if (classif.method == "nsc") {

    library(pamr); 

    # Nearest shrunken centroids

    train.data    <- list(x = emat[, train.samples], y = factor(cl[train.samples]))
    test.data     <- list(x = emat[, test.samples],  y = factor(cl[test.samples]))

    # PAMR and CV
    res  <- pamr.train(train.data)
    mycv <- pamr.cv(res, train.data, nfold = 10)

    # Classify new samples
    pred <- pamr.predict(res, emat[, test.samples], threshold = mycv$threshold[which.min(mycv$error)], type = "class")

    # Prediction accuracy
    acc <- mean(pred == cl[test.samples])

    accs[cvn, method] <- acc

  } else if (classif.method == "random.forest") {

     # Random forest
     library(randomForest)

     # Training data
     train.data <- data.frame(cbind(t(emat[gene.subset, train.samples])))
     train.data$classInfo <- factor(cl[train.samples])

     # Test data
     test.data <- t(emat[gene.subset, test.samples])
     colnames(test.data) <- paste("X", colnames(test.data), sep = "")

     # Fit the model
     rf <- randomForest(classInfo ~ ., data = train.data, importance = FALSE, proximity = FALSE, ntree = ntree)

     # Test prediction on test data
     pred <- predict(rf, test.data)		  
     acc <- mean(pred == cl[test.samples])

     accs[cvn, method] <- acc

  }

  print(paste(cvn, method, acc))

}

}


# Testing if row method is less accurate than col method
pvals <- matrix(NA, nrow = ncol(accs), ncol = ncol(accs))
rownames(pvals) <- colnames(pvals) <- colnames(accs)
for (method1 in colnames(accs)) {
  for (method2 in colnames(accs)) {
    pvals[method1, method2] <- wilcox.test(accs[, method1], accs[, method2], "less", paired = TRUE)$p.value
  }
}

write.table(round(pvals,5), sep = "\t", file = "Classification-comparisons-ntree500fs01.tab")


save(pvals, accs, class.minsize, fs, ntree, file = "accs.minsize5.RData")



library(ggplot2)
library(reshape)
o <- order(apply(accs, 2, mean))
accs <- accs[,o]
df <- data.frame(accs); df$fold <- 1:N
df <- melt(df, id = c("fold"))
df$method <- df$variable
p <- ggplot(df) + aes(x = method, y = value, group = method, fill = method) + geom_boxplot()

pdf("Classification-comparisons-ntree500fs01.pdf")
print(p)
dev.off()

