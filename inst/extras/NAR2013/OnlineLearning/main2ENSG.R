print("Sample classes")
load("data/classInfo.RData") # classInfo

emat.files <- c(rpa.online = "data/emat.rpa.online.storage.ENSG.RData",  
	        #rpa.online.berbg = "data/emat.rpa.online.berbg.ENSG.RData",  
	        #rpa.online.bersd = "data/emat.rpa.online.bersd.ENSG.RData",  
	        #rpa.online.bermean = "data/emat.rpa.online.bermean.ENSG.RData",  
	        #rpa.online.lm = "data/emat.rpa.online.lm.ENSG.RData",  
	        rma = "data/Lukk-RMA-HGU133A_Hs_ENSG.RData", 
		frma = "data/Lukk-fRMA-HGU133A_Hs_ENSG.RData", 
		frmab = "data/Lukk-fRMA.batch-HGU133A_Hs_ENSG.RData", 
		mas = "data/Lukk-MAS5-HGU133A_Hs_ENSG.RData")


rs <- names(classInfo)
classInfo <- classInfo[rs]

# Remove smallest classes
tab <- table(classInfo)

accs.list <- list()

#minsizes <- c(1, 2, 3)
minsizes <- 2 # singletons excluded

#set.seed(342)
set.seed(3422)
N <- 10
ntree <- 500
#fs <- 0.1
classif.method <- "random.forest"

# Random forest
library(randomForest)
source("~/scripts/R/Affy.R")

for (class.minsize in minsizes) {

  selected.classes <- names(which(tab >= class.minsize) )
  s <- names(classInfo)[which(classInfo %in% selected.classes)]
  cl <- classInfo[s]

  # Pick one of the data matrices
  load(emat.files[["rma"]])
  emat <- emat.rma.ENSG[, s]

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
    all.sets <- removeAFFX(rownames(emat.rma.ENSG))
    #gene.subset <- sample(all.sets, round(fs*length(all.sets)))
    gene.subset <- sample(all.sets, 1000)

    for (method in names(emat.files)) {

      print(paste(class.minsize, cvn, method))

      # Pick data for the selected method
      load(emat.files[[method]])  
      if (method == "rpa") { emat <- emat.rpa.ENSG[, s]; rm(emat.rpa.ENSG)}    # version used in the other experiments
      if (method == "rpa.online") { emat <- emat.rpa.online.storage.ENSG[, s]; rm(emat.rpa.online.storage.ENSG)} # latest version
      if (method == "rpa.online.berbg") { emat <- emat.rpa.online.berbg.ENSG[, s]; rm(emat.rpa.online.berbg.ENSG)} # latest version
      if (method == "rpa.online.bersd") { emat <- emat.rpa.online.bersd.ENSG[, s]; rm(emat.rpa.online.bersd.ENSG)} # latest version
      if (method == "rpa.online.bermean") { emat <- emat.rpa.online.bermean.ENSG[, s]; rm(emat.rpa.online.bermean.ENSG)} # latest version
      if (method == "rpa.online.lm") { emat <- emat.rpa.online.lm.ENSG[, s]; rm(emat.rpa.online.lm.ENSG)} # latest version
      if (method == "rma") { emat <- emat.rma.ENSG[, s]; rm(emat.rma.ENSG)}
      if (method == "frma") { emat <- emat.frma.ENSG[, s]; rm(emat.frma.ENSG)}
      if (method == "frmab") { emat <- emat.frma.batch.ENSG[, s]; rm(emat.frma.batch.ENSG)}
      if (method == "mas") { emat <- emat.mas5.ENSG[, s]; rm(emat.mas5.ENSG)}

      # Training data
      train.data <- data.frame(cbind(t(emat[gene.subset, train.samples])))
      train.data$classInfo <- factor(cl[train.samples])

      # Test data
      test.data <- t(emat[gene.subset, test.samples])
      #colnames(test.data) <- paste("X", colnames(test.data), sep = "")

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


save(minsizes, accs.list, ntree, file = "Classification-comparisons-ENSG.RData")

