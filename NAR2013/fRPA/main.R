# Get the data (E-GEOD-5258)
# Transcription profiling by array of human cell lines after treatment with various small molecules
# The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. Justin Lamb, Emily D Crawford, David Peck, Joshua W Modell, Irene C Blat, Matthew J Wrobel, Jim Lerner, Jean-Philippe Brunet, Aravind Subramanian, Kenneth N Ross, Michael Reich, Haley Hieronymus, Guo Wei, Scott A Armstrong, Stephen J Haggarty, Paul A Clemons, Ru Wei, Steven A Carr, Eric S Lander, Todd R Golub. , Europe PMC 17008526

# Load the data
#system("wget http://www.ebi.ac.uk/arrayexpress/files/E-GEOD-5258/E-GEOD-5258.eSet.r")
load("E-GEOD-5258.eSet.r")
library(affy)
abatch <- study[[2]]

# -----------------------------------------------------------

# Expression matrices
source("get.esets.R")

# -----------------------------------------------------------

# Get classifications
#classInfo <- phenoData(study[[2]])@data[["Factor.Value..CellLine."]]
classInfo <- phenoData(study[[2]])@data[["Factor.Value..Compound."]]
#classInfo <- phenoData(study[[2]])@data[["Factor.Value..Vehicle."]]
names(classInfo) <- rownames(phenoData(study[[2]])@data)

class.minsize <- 2 # minimal classes excluded
set.seed(3421)
N <- 5
ntree <- 200
classif.method <- "random.forest"

# Random forest
library(randomForest)
source("~/scripts/R/Affy.R")

# Remove smallest classes
tab <- table(classInfo)
selected.classes <- names(which(tab >= class.minsize) )
s <- names(classInfo)[which(classInfo %in% selected.classes)]
cl <- classInfo[s]
all.sets <- removeAFFX(rownames(esets[["RMA"]]))

# Initialize randomforest table
accs <- matrix(NA, nrow = N, ncol = length(esets))
colnames(accs) <- names(esets)

# Split data in training and test sets N-fold
folds <- list()
for (cvn in 1:N) {
  n <- floor(length(cl)/N)
  if (cvn == N) {n <- length(s)}
  inds <- sample(s, n)
  folds[[cvn]] <- inds
  s <- setdiff(s, inds)
}


set.seed(553)
for (cvn in 1:N) {

    # Split training/test data (small training data, large test data)
    train.samples <- unlist(folds[cvn])
    test.samples  <- unlist(folds[-cvn])

    # Select random subset of genes to speed up
    gene.subset <- sample(all.sets, 1000)

    for (method in names(esets)) {

      print(paste(class.minsize, cvn, method))
      emat <- esets[[method]]
      colnames(emat) <- gsub(".gz", "", colnames(emat))
      colnames(emat) <- gsub(".CEL.gz", "", colnames(emat))
      colnames(emat) <- sapply(strsplit(colnames(emat), "_"), function (x) {x[[1]]})

      # Training data
      train.data <- data.frame(cbind(t(emat[gene.subset, train.samples])))
      train.data$classInfo <- factor(cl[sapply(strsplit(train.samples, "_"), function (x) {x[[1]]})])

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

boxplot(accs[, order(apply(accs, 2, median))])







  df <- data.frame(accs); 
  df$fold <- 1:nrow(accs)
  library(reshape)
  df <- melt(df, id = c("fold"))
  df$method <- df$variable
  df$method <- factor(df$method, levels = c("MAS5", "RMA", "fRMA", "RPA", "fRMA.batch"))


  # Testing if row method is less accurate than col method
  pvals <- matrix(NA, nrow = ncol(accs), ncol = ncol(accs))
  rownames(pvals) <- colnames(pvals) <- colnames(accs)
  for (method1 in colnames(accs)) {
    for (method2 in colnames(accs)) {
      pvals[method1, method2] <- wilcox.test(accs[, method1], accs[, method2], "less", paired = TRUE)$p.value
    } 
  }    
  print(round(pvals, 5))
  #write.table(round(pvals,5), sep = "\t", file = paste("Classification-comparisons.", maintext, k, ".tab", sep = ""))

  library(ggplot2)
  p <- ggplot(df) + aes(x = method, y = value) + theme_set(theme_bw(20)) + ylab("Classification accuracy") + xlab("") + ggtitle("Random forest comparison") + opts(axis.text.x = theme_text(angle = 25, size = 20)) + opts(legend.position = "none") + geom_boxplot() 
  #p <- p + scale_y_continuous(limits = lims))
  # opts(axis.text.x = theme_text(vjust = 0.5)

  print(p)


