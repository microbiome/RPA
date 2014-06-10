# Preprocess data set with various methods
esets <- list()

print("fRPA")
library(RPA)
# hgu133a.rpa.priors
load("~/Rpackages/RPA/github/RPA/OnlineLearning/HGU133A-RPA-priors.RData")

esets$fRPA <- exprs(frpa(abatch, probe.parameters = hgu133a.rpa.priors))
colnames(esets$fRPA) <- gsub(".CEL", "", colnames(esets$fRPA))
colnames(esets$fRPA) <- gsub(".gz", "", colnames(esets$fRPA))

print("RPA")
library(RPA)
esets$RPA <- exprs(rpa(abatch))
colnames(esets$RPA) <- gsub(".CEL", "", colnames(esets$RPA))

print("fRMA")
library(frma)
esets$fRMA <- exprs(frma(abatch))

print("RMA")
esets$RMA <- exprs(rma(abatch))

#print("MAS")
esets$MAS <- log2(exprs(mas5(abatch)))


# Test
#hgu133afrmavecs
eset.frma <- frma(abatch, input.vecs = hgu133afrmavecs)
