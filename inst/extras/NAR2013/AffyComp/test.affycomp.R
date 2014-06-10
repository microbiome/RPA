#############################################

# Preprocess spike-in data sets with RPA

require(affy)
library(RPA)
#source("http://bioconductor.org/biocLite.R"); biocLite("hgu133atagcdf")
library("hgu133atagcdf")

# get affybatch
data.path <- "hgu133spikein"
CELlist <- list.celfiles(data.path, full.names = TRUE)
abatch <- ReadAffy( filenames = CELlist )

# Get set.inds for HGU133Atag
pN <- probeNames(abatch)
set.inds <- split(1:length(pN), pN) # pNList

res <- rpa.complete(abatch)
#ptab <- probetable(res$probe.parameters)

# hgu133a.rpa.priors
load("~/Rpackages/RPA/github/RPA/OnlineLearning/HGU133A-RPA-priors.RData")
#hgu133a.rpa.priors
hgu133a.priors.list <- probe.parameters.tolist(hgu133a.rpa.priors)

# Initialise fRPA parameters
probe.parameters2 <- res$probe.parameters

common.sets <- intersect(names(probe.parameters2$tau2), names(hgu133a.priors.list$tau2)) 
for (set in common.sets) {
  probe.parameters2$tau2[[set]] <- hgu133a.priors.list$tau2[[set]]
  probe.parameters2$affinity[[set]] <- hgu133a.priors.list$affinity[[set]]
} 

probe.parameters2$quantile.basis[na.omit(match(names(unlist(probe.parameters2$tau2)), names(unlist(hgu133a.priors.list$tau2))))] <- hgu133a.priors.list$quantile.basis
probe.parameters2$alpha <- NULL
probe.parameters2$betas <- NULL

# Original probe parameters table
tabo <- probetable(res$probe.parameters)
# fRPA probe parameters table
tabf <- probetable(probe.parameters2)

rs <- sample(nrow(tabo), 1e3)
par(mfrow = c(2,2))
for (varname in c("quantile.basis", "tau2", "affinity")) {
  plot(tabo[rs, varname], tabf[rs, varname], pch = ".", xlab = "Spikein", ylab = "Lukk", main = varname); abline(0,1,lty = 2)
}


emat <- exprs(res$eset)
topsets <- rownames(emat)[rev(order(apply(emat, 1, sd)))[1:50]]
par(mfrow = c(2,2))
for (varname in c("quantile.basis", "tau2", "affinity")) {
  plot(subset(tabo, probeset %in% topsets)[, varname], subset(tabf, probeset %in% topsets)[, varname], pch = ".", xlab = "Spikein", ylab = "Lukk", main = varname); abline(0,1,lty = 2)
}


varname <- "tau2"
k <- k+1; set <- topsets[[k]]; plot(sqrt(subset(tabo, probeset %in% set)[, varname]), sqrt(subset(tabf, probeset %in% set)[, varname]), xlab = "Original", ylab = "Lukk", main = varname); abline(0,1)

# -------------------------------------------

x <- exprs(res$expressionSet)
colnames(x) <- sapply(colnames(x), function (x){strsplit(x, "\\.")[[1]][[1]]})
write.table(data.frame(2^x,check.names=FALSE),file=paste(data.path, ".csv", sep = ""),sep=",",col.names=NA,quote=FALSE)


