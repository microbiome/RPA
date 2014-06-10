# With respect to the fRMA, if you create a vector of CEL file names, 
# 'CELlist', this is the code you could use to preprocess the data:

##################################################

# necessary packages

source("get.hgu133atag.priors.R")

# -----------------

eset <- exprs(rpa(cel.files = CELlist, probe.parameters = hgu133atag.rpa.priors.list, summarize.with.affinities = TRUE))
 nam <- "fRPA"
 x <- eset
 colnames(x) <- sapply(colnames(x), function (x){strsplit(x, "\\.")[[1]][[1]]})
 write.table(data.frame(2^x,check.names=FALSE), file = paste(nam, ".csv", sep = ""),sep=",",col.names=NA,quote=FALSE)

# -> not getting fRPA working with spikein data. perhaps since the
#     extra probesets are informative or since spikein data sets in
#     general are very different from real arrays?

# 17 AFFX control sets are different
# hgu133afrmavecs
# hgu133atagfrmavecs
