# With respect to the fRMA, if you create a vector of CEL file names, 
# 'CELlist', this is the code you could use to preprocess the data:

##################################################

# necessary packages

require(affy)
require(frma)
require(hgu133afrmavecs)
data(hgu133afrmavecs)
library(RPA)

data.path <- "hgu133spikein"

# run fRMA, sample by sample
CELlist <- list.celfiles(data.path, full.names = TRUE)
myData <- ReadAffy( filenames = CELlist )
#myData@cdfName <- "hgu133atag"

esets <- list()
esets$RPA <- exprs(rpa(cel.files = CELlist))
esets$RMA <- exprs(justRMA(filenames = CELlist))

# Get set.inds for HGU133Atag
#myData@cdfName <- "hgu133a"
pN <- probeNames(myData)
set.inds <- split(1:length(pN), pN) # pNList

# Get HGU133A abatch
# Su et al 2004 data
hgu133a.cels <- sample(list.celfiles("/home/BACKUPS/Transcend3/data/Affy/Su2004/CEL/hgu133a", full.names = TRUE), 5)
abatch.hgu133a <- ReadAffy(filenames = hgu133a.cels)
pN.hgu133a <- probeNames(abatch.hgu133a)
set.inds.hgu133a <- split(1:length(pN.hgu133a), pN.hgu133a) # pNList

# 17 AFFX control sets are different
# hgu133afrmavecs
# hgu133atagfrmavecs

# -----------------------------------------------------

# Initialize fRMA and RPA vector lists
hgu133atagfrmavecs <- list()
# Initialize probe vector
pvec <- lapply(set.inds, function (x) {rep(NA, length(x))}); names(pvec) <- names(set.inds)
# Order sets in the same order as in hgu133a, then the rest
pvec <- pvec[c(names(set.inds.hgu133a), setdiff(names(pvec), names(set.inds.hgu133a)))]
pvec.init <- unlist(pvec)
for (varname in names(hgu133afrmavecs)) {

  # Copy the contents from hgu133a
  pvec <- pvec.init
  pvec[1:length(hgu133afrmavecs[[varname]])] <- hgu133afrmavecs[[varname]]
  hgu133atagfrmavecs[[varname]] <- pvec
  hgu133atagfrmavecs[[varname]][is.na(hgu133atagfrmavecs[[varname]])] <- median(na.omit(hgu133atagfrmavecs[[varname]]))

}

# -----------------------------------------------------

# hgu133a.rpa.priors
load("~/Rpackages/RPA/github/RPA/OnlineLearning/HGU133A-RPA-priors.RData")

# Initialize probe vector
#pp <- probe.parameters.tolist(hgu133a.rpa.priors)

# Initialize probe vector
pvec <- lapply(set.inds, function (x) {rep(NA, length(x))}); 
names(pvec) <- names(set.inds)
pvec <- pvec[c(names(set.inds.hgu133a), setdiff(names(pvec), names(set.inds.hgu133a)))]
ptab <- list(tau2 = unlist(pvec), probeset = unlist(sapply(1:length(pvec), function (k) { rep(names(pvec)[[k]], length(pvec[[k]])) })), probe.index = unlist(sapply(1:length(pvec), function (k) { 1:length(pvec[[k]]) })))
ptab$tau2[1:length(hgu133a.rpa.priors$tau2)] <- hgu133a.rpa.priors$tau2 
ptab$quantile.basis <- rep(NA, length(ptab$tau2))
ptab$quantile.basis[1:length(ptab$quantile.basis)] <- hgu133a.rpa.priors$quantile.basis 
ptab$quantile.basis[is.na(ptab$quantile.basis)] <- median(na.omit(ptab$quantile.basis))
ptab <- data.frame(ptab)
#esets$fRPA <- rpa(cel.files = CELlist, probe.parameters = ptab)
pvecs <- probe.parameters.tolist(ptab)
esets$fRPA <- exprs(rpa(cel.files = CELlist, probe.parameters = pvecs))
# -> not getting fRPA working with spikein data. perhaps since the extra probesets are informative?

esets$fRMA <- exprs(frma(myData, summarize = "robust_weighted_average", input.vecs = hgu133atagfrmavecs))


#rs <- sample(1e5, 1e3); plot(exprs(eset.rma)[rs], exprs(eset.fRMA)[rs], pch = "."); abline(0,1)
#rs <- sample(1e5, 1e3); plot(esets$RPA[rs], esets$fRPA[rs], pch = "."); abline(0,1)

##################################################

for (nam in names(esets)) {
 x <- esets[[nam]]
 colnames(x) <- sapply(colnames(x), function (x){strsplit(x, "\\.")[[1]][[1]]})
 write.table(data.frame(2^x,check.names=FALSE), file = paste(nam, ".csv", sep = ""),sep=",",col.names=NA,quote=FALSE)

}



# -> missing some of the spikes ("Tag") -> cannot be tested
#library("pd.hg.u133a.tag")
#myData.fRMA <- frma(myData, summarize = "robust_weighted_average")
