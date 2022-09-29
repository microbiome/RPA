# Create ber versions first with batch.effects.R

print("Sample classes")
load("data/classInfo.RData") # classInfo

emat.files <- c(rpa.online = "data/emat.rpa.online.storage.RData",  
	        #rpa.online2 = "~/Rpackages/RPA/github/RPA/OnlineLearning/LukkOutput20121104.RData", # improved affinity thing
	        #rpa.online.berbg = "data/emat.rpa.online.berbg.RData",  
	        #rpa.online.bersd = "data/emat.rpa.online.bersd.RData",  
	        #rpa.online.bermean = "data/emat.rpa.online.bermean.RData",  
	        #rpa.online.lm = "data/emat.rpa.online.lm.RData",  
	        rma = "data/Lukk-RMA.RData", 
		frma = "data/Lukk-fRMA.RData", 
		frmab = "data/Lukk-fRMA.batch.RData",
		mas = "data/Lukk-MAS5.RData"
		)

rs <- names(classInfo)
classInfo <- classInfo[rs]

# -----------------------------------------------------------------------------------------------

accs.list <- list()

# minsizes <- c(1, 2, 3)
minsizes <- 2 # singletons excluded

set.seed(3421)
N <- 10
ntree <- 500
#fs <- 0.1
classif.method <- "random.forest"

# Random forest
library(randomForest)
source("~/scripts/R/Affy.R")

for (class.minsize in minsizes) {

  # Remove smallest classes
  tab <- table(classInfo)
  selected.classes <- names(which(tab >= class.minsize) )
  s <- names(classInfo)[which(classInfo %in% selected.classes)]
  cl <- classInfo[s]

  # Pick one of the data matrices
  load(emat.files[["rma"]])
  emat <- emat.rma[, s]

  accs <- matrix(NA, nrow = N, ncol = length(emat.files))
  colnames(accs) <- names(emat.files)

  # Split data in training and test sets N-fold
  folds <- list()
  samples <- colnames(emat)
  for (cvn in 1:N) {
    n <- floor(ncol(emat)/N)
    if (cvn == N) {n <- length(samples)}
    inds <- sample(samples, n)
    folds[[cvn]] <- inds
    samples <- setdiff(samples, inds)
  }

  for (cvn in 1:N) {

    gc()

    # Split training/test data
    train.samples <- unlist(folds[-cvn])
    test.samples  <- unlist(folds[cvn])

    load(emat.files[["rma"]])
    # Select random subset of genes to speed up
    all.sets <- removeAFFX(rownames(emat.rma))
    #gene.subset <- sample(all.sets, round(fs*length(all.sets)))
    gene.subset <- sample(all.sets, 1000)

    for (method in names(emat.files)) {

      print(paste(class.minsize, cvn, method))

      # Pick data for the selected method
      load(emat.files[[method]])
      if (method == "rpa") { emat <- emat.rpa[, s]; rm(emat.rpa)}    # version used in the other experiments
      if (method == "rpa.online") { emat <- emat.rpa.online.storage[, s]; rm(emat.rpa.online.storage)} # latest version
      if (method == "rpa.online2") { emat <- emat[, s]; rm(emat.rpa.online.storage)} # new affinity version
      if (method == "rpa.online.berbg") { emat <- emat.rpa.online.berbg[, s]; rm(emat.rpa.online.berbg)} # latest version
      if (method == "rpa.online.bersd") { emat <- emat.rpa.online.bersd[, s]; rm(emat.rpa.online.bersd)} # latest version
      if (method == "rpa.online.bermean") { emat <- emat.rpa.online.bermean[, s]; rm(emat.rpa.online.bermean)} # latest version
      if (method == "rpa.online.lm") { emat <- emat.rpa.online.lm[, s]; rm(emat.rpa.online.lm)} 
      if (method == "rma") { emat <- emat.rma[, s]; rm(emat.rma)}
      if (method == "frma") { emat <- emat.frma[, s]; rm(emat.frma)}
      if (method == "frmab") { emat <- emat.frma.batch[, s]; rm(emat.frma.batch)}
      if (method == "mas") { emat <- emat.MAS5[, s]; rm(emat.MAS5)}

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

      print(acc)

    }

  } 

  accs.list[[as.character(class.minsize)]] <- accs

}

save(minsizes, accs.list, ntree, file = "Classification-comparisons.RData")


