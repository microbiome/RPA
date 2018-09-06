#if (!requireNamespace("BiocManager", quietly=TRUE))
    #install.packages("BiocManager")
#BiocManager::install("hgu133plus2hsensgcdf")
#library("hgu133plus2hsensgcdf")
#cel.files = cels; mc.cores = 4; batch.size = 4; batch.file.id = "batch"; 
#cdf = cdf; cel.path = NULL; sets = NULL; bg.method = "rma"; 
#priors = list(alpha = 2, beta = 1); epsilon = 1e-2; cind = 1; verbose = TRUE; 
#shuffle = TRUE; batches = NULL; 
#batch.list.file.id = NULL; quantile.basis = NULL; quantile.file.id = NULL; 
#hyper.parameters = NULL;  hyperparameter.batch.id = NULL 

set.seed(11122) 
require(RPA)

#fs <- list.files("../RPA/R", full.names = TRUE, pattern = ".R$")
#for (f in fs) {source(f)}

#cdf <- NULL
cdf <- "HGU133Plus2_Hs_ENSG"

cels <- list.celfiles("/share/mi/data/GSE3526/CEL/", full.names = T)[1:50]

st <- Sys.time()
#eset <-  rpa.online(cel.files = cels, mc.cores = 4, batch.size = 4, batch.file.id = "batch", cdf = cdf)
eset <- rpa.online(cel.files = cels, mc.cores = 4, batch.size = 25, cdf = cdf, save.batches = TRUE)
et <- Sys.time()

# OUTOA: batch-filujen kanssa saa paremmat korrelaatiot RPA vs RPA-online? 

print(et-st)

emat <- exprs(eset)
colnames(emat) <- unname(sapply(colnames(emat), function(s){unlist(strsplit(s, "/"))[[8]]}))

##########################


check <- TRUE
if (check) {

  # Compare to Original RPA and RMA

  # Compare
  abatch.tmp <- ReadAffy(filenames = cels)
  if (!is.null(cdf)) { abatch.tmp@cdfName <- cdf }
  
  # rpa <- RPA.pointestimate(abatch.tmp, cdf = cdf)

  eset.rma <- exprs(rma(abatch.tmp))
  eset.rpa <- exprs(rpa(abatch.tmp, cdf = cdf))
  s <- colnames(eset.rma)
  save(s, eset.rma, eset.rpa, file = "tmp.RData")
}
load("tmp.RData")


#########################

  cors <- matrix(NA, nrow = nrow(eset.rma), ncol = 3)
  rownames(cors) <- rownames(eset.rma)
  colnames(cors) <- c("RMA.RPA", "RMA.RPAonline", "RPA.RPAonline")
  for (set in sample(rownames(eset.rma), 1000)) {
    print(set)
    cors[set,] <- c(cor(eset.rma[set,s], eset.rpa[set,s]),
    cor(eset.rma[set,s], emat[set,s]),
    cor(eset.rpa[set,s], emat[set,s]))
  }

par(mfrow = c(3,1))
for (i in 1:ncol(cors)) {
  hist(cors[,i], 100, main = colnames(cors)[[i]], xlim = c(-1,1))
}

k <- 1
df <- data.frame(cor(cbind(rma = eset.rma[,s[[k]]], rpa = eset.rpa[,s[[k]]], rpao = emat[,s[[k]]])))
print(df)

#print(et-st)






