#############################################

# With RPA_1.7.33 and R-2.12.2

#############################################

# Preprocess spike-in data sets with RPA

library(RPA)


#for (data.path in c("hgu133spikein", "spikein")) {
for (data.path in c("hgu133spikein")) {

  abatch <- ReadAffy(celfile.path = data.path)

  x <- exprs(rpa(abatch))
  #eset.rma <- rma(abatch)

  colnames(x) <- sapply(colnames(x), function (x){strsplit(x, "\\.")[[1]][[1]]})

  write.table(data.frame(2^x,check.names=FALSE),file=paste(data.path, ".csv", sep = ""),sep=",",col.names=NA,quote=FALSE)

}

###############################################
