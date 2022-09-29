library(affy)

# Su et al 2004 data
cels <- list.celfiles("/home/BACKUPS/Transcend3/data/Affy/Su2004/CEL/hgu133a", full.names = TRUE)
abatch <- ReadAffy(filenames = cels)
library(gdata)
ann <- read.xls("/home/BACKUPS/Transcend3/data/Affy/Su2004/CEL/hgu133a/su2004tissues.xls")

file.ids <- lapply(as.character(ann$X3A.filenames), function (x) strsplit(x, ",")[[1]])
names(file.ids) <- as.character(ann$Tissue)
file.ids <- unlist(file.ids)
ids <- sapply(names(file.ids), function (x) {substr(x, 1, nchar(x)-1)})
names(ids) <- file.ids
classInfo <- ids
# vrt. colnames(exprs(abatch))

# -----------------------------------------

source("get.esets.R")

classInfo <- gsub("isletcell","islet",substr(gsub("-", "", gsub(" ", "", gsub("_", "", tolower(gsub(".CEL.gz", "", sapply(strsplit(colnames(esets$RMA), "_"), function (x) {paste(x[-1], collapse = "_")})))))), 1, 12))
names(classInfo) <- gsub(".CEL.gz", "", colnames(esets$RMA))

training.samples <- gsub(".CEL.gz", "", names(classInfo)[match(unique(classInfo), classInfo)])
test.samples <- gsub(".CEL.gz", "", setdiff(names(classInfo), training.samples))

training.classes <- classInfo[training.samples]
test.classes <- classInfo[test.samples]


# -----------------------------------------------------------

uniq.classes <- unique(classInfo)

eset <- esets$RMA

set.seed(3533)
N <- 1000
accmat <- matrix(NA, nrow = N, ncol = 5)
colnames(accmat) <- names(esets)
for (n in 1:N) {
  print(n)
  accs <- c()
  sample.cors <- matrix(NA, ncol = length(esets), nrow = length(uniq.classes))
  colnames(sample.cors) <- names(esets)
  rownames(sample.cors) <- uniq.classes
  sample.cors.rand <- sample.cors

  for (method in names(esets)) {

  eset <- esets[[method]]
  colnames(eset) <- gsub(".CEL.gz", "", colnames(eset))
  colnames(eset) <- gsub(".gz", "", colnames(eset))

  #eset.train <- eset[, training.samples]
  #colnames(eset.train) <- classInfo[training.samples]
  #eset.test <- eset[, test.samples]
  #colnames(eset.test) <- classInfo[test.samples]

  # bootstrap sample of probesets
  #s <- sample(nrow(eset), replace = TRUE)
  # Jackknife 
  s <- sample(nrow(eset), floor(nrow(eset)*0.2))
  eset.train <- eset[s, training.samples]
  colnames(eset.train) <- classInfo[training.samples]
  eset.test <- eset[s, test.samples]
  colnames(eset.test) <- classInfo[test.samples]

  eset.test <- eset.test[, colnames(eset.train)]
  eset.test.rand <- eset.test
  colnames(eset.test.rand) <- sample(colnames(eset.test))

  for (cl in uniq.classes) { sample.cors[cl, method] <- cor(eset.train[, cl], eset.test[, cl], method = "spearman")  }
  for (cl in uniq.classes) { sample.cors.rand[cl, method] <- cor(eset.train[, cl], eset.test.rand[, cl], method = "spearman")  }

  cc <- cor(eset[s, test.samples], eset[s, training.samples], method = "spearman")
  #cc <- cor(eset[s, test.samples], eset[s, training.samples], method = "pearson")
  rownames(cc) <- classInfo[rownames(cc)]
  colnames(cc) <- classInfo[colnames(cc)]

  cc <- cc[uniq.classes, uniq.classes]
  # Check that the correct class has the maximal similarity
  acc <- c(); for (k in 1:nrow(cc)) {acc[[k]] <- max(cc[k, ]) == cc[k,k]}
  accs[[method]] <- mean(acc)

  }

  accmat[n, names(accs)] <- accs

}



o <- order(apply(accmat[, setdiff(colnames(accmat), "MAS")], 2, mean))
boxplot(accmat[, o], las = 2)

#library(ggplot2); library(reshape); df <- melt(accmat); p <- ggplot(df, aes(x = X2, y = value)) + geom_boxplot(); print(p)


# Row method greater than col method?
pvals.accs <- matrix(NA, nrow = ncol(accmat), ncol = ncol(accmat))
rownames(pvals.accs) <- colnames(pvals.accs) <- colnames(accmat)
for (method1 in colnames(accmat)) {
  for (method2 in setdiff(colnames(accmat), method1)) {    
    pvals.accs[method1, method2] <- wilcox.test(accmat[, method1], accmat[, method2], alternative = "greater", paired = TRUE)$p.value
  }
}

print(pvals.accs)

par(mfrow = c(4,1))
for (nam in sort(colnames(accmat))) {
  hist(accmat[,nam], 100, xlim = range(accmat), main = nam)
}

# ----------------------------


#sort(apply(accmat, 2, median))
#wilcox.test(accmat[, "fRPA"], accmat[, "fRMA"], alternative = "greater", paired = TRUE)

#par(mfrow = c(1, 2))
#o <- order(apply(sample.cors, 2, median))
#boxplot(sample.cors[, o], las = 2)
#boxplot(sample.cors.rand[, o], las = 2)

# Row method greater than col method?
pvals <- matrix(NA, nrow = ncol(sample.cors), ncol = ncol(sample.cors))
rownames(pvals) <- colnames(pvals) <- colnames(sample.cors)
for (method1 in colnames(sample.cors)) {
  for (method2 in colnames(sample.cors)) {
    pvals[method1, method2] <- wilcox.test(sample.cors[, method1], sample.cors[, method2], alternative = "greater", paired = TRUE)$p.value
  }
}

# Row method greater than col method?
pvals.rand <- matrix(NA, nrow = ncol(sample.cors), ncol = ncol(sample.cors))
rownames(pvals.rand) <- colnames(pvals.rand) <- colnames(sample.cors)
for (method1 in colnames(sample.cors)) {
  for (method2 in colnames(sample.cors)) {
    pvals.rand[method1, method2] <- wilcox.test(sample.cors.rand[, method1], sample.cors.rand[, method2], alternative = "greater")$p.value
  }
}



