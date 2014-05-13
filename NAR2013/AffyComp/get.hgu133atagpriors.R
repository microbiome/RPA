



# Get HGU133A abatch
# Su et al 2004 data
library(affy)
#hgu133a.cels <- sample(list.files("/home/BACKUPS/Transcend3/data/Affy/Su2004/CEL/hgu133a", full.names = TRUE, pattern = "CEL.gz"), 5)

abatch.hgu133a <- ReadAffy(filenames = hgu133a.cels)
pN.hgu133a <- probeNames(abatch.hgu133a)
set.inds.hgu133a <- split(1:length(pN.hgu133a), pN.hgu133a) # pNList

# hgu133a.rpa.priors
load("~/Rpackages/RPA/github/RPA/OnlineLearning/HGU133A-RPA-priors.RData")

# Initialize probe vector
#pp <- probe.parameters.tolist(hgu133a.rpa.priors)

# Initialize probe vector
hgu133atag.rpa.priors <- list()
hgu133atag.rpa.priors$quantile.basis <- hgu133a.rpa.priors$quantile.basis[match(unlist(set.inds), unlist(set.inds.hgu133a))]
# Replace missing with median
hgu133atag.rpa.priors$quantile.basis[is.na(hgu133atag.rpa.priors$quantile.basis)] <- median(na.omit(hgu133atag.rpa.priors$quantile.basis))

hgu133atag.rpa.priors$tau2 <- hgu133a.rpa.priors$tau2[match(unlist(set.inds), unlist(set.inds.hgu133a))]
# Replace missing with uninformative prior
hgu133atag.rpa.priors$tau2[is.na(hgu133atag.rpa.priors$tau2)] <- 1

hgu133atag.rpa.priors$affinity <- hgu133a.rpa.priors$affinity[match(unlist(set.inds), unlist(set.inds.hgu133a))]
# Replace missing with uninformative prior
hgu133atag.rpa.priors$affinity[is.na(hgu133atag.rpa.priors$affinity)] <- 1

hgu133atag.rpa.priors$probeset <- unlist(sapply(names(set.inds), function (set) {rep(set, length(set.inds[[set]]))}))
hgu133atag.rpa.priors$probe.index <- unlist(sapply(names(set.inds), function (set) {1:length(set.inds[[set]])}))
hgu133atag.rpa.priors <- data.frame(hgu133atag.rpa.priors)
hgu133atag.rpa.priors.list <- probe.parameters.tolist(hgu133atag.rpa.priors)

#library("pd.hg.u133a.tag")
#data("pd.hg.u133a.tag")
#myData@cdfName <- "hgu133atag"
#myData@cdfName <- "hgu133a"
#eset <- exprs(rpa(cel.files = CELlist))
