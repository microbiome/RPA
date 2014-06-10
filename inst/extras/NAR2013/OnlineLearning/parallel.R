#install.packages("Rpackages/RPA_1.3.3.tar.gz")                                        

require(affy)
require(RPA)
require(affyPara)
require(hgu133plus2cdf)
require(hgu133plus2hsentrezgcdf)

# List CEL files
tumor.path <- "/wrk/lmlahti/GSE2109"
cels.t <- list.celfiles(tumor.path, full.names = TRUE)
normal.path <- "/wrk/lmlahti/GSE3526"
cels.n <- list.celfiles(normal.path, full.names = TRUE)
cels <- c(cels.t, cels.n)
#cels <- cels.n

# Start a 'cluster' (multiples of 12 in Triton)
cl <- makeCluster(2 * 12, type = "SOCK")

# NOTE: this assumes that all CEL files are available for all cluster nodes
# background correct
abatch <- bgCorrectPara(cels, method = "rma", verbose = TRUE)
# normalize
abatch <- normalizeAffyBatchQuantilesPara(abatch, type = "pmonly", verbose = TRUE)

# Preprocess
#eset <- rpa(abatch, cdfname = "HGU133Plus2_Hs_ENTREZG", verbose = TRUE)
eset <- rpa(abatch, cdfname = "", verbose = TRUE)
save(eset, file="/wrk/lmlahti/GSE2109andGSE3526-RPA-altcdf.Rdata")

stopCluster()

print(sessionInfo())

