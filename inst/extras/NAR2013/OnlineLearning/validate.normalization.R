
require(RPA)
require(multicore)
fs <- list.files("../RPA/R", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}

set.seed(11122)

 	        cel.path = NULL
cel.files = list.celfiles("/share/mi/data/GSE3526/CEL/", full.names = T)[1:12]
                    sets = NULL
                  priors = list(alpha = 2, beta = 1)
                 epsilon = 1e-2 
                    cind = 1
                 verbose = FALSE
               bg.method = "rma"    
    normalization.method = "quantiles"
                     cdf = NULL 
              batch.size = 4
	  quantile.basis = NULL 
           batch.file.id = "batch.id"

# Ordinary quantile-normalized data
abatch <- ReadAffy(filenames = cel.files)
abatch.bg <- bg.correct(abatch, bg.method, destructive = TRUE)
ab.q2 <- do.call(affy:::normalize, c(alist(abatch.bg, "quantiles"), normalize.param = list()))
qmat2 <- log2(pm(ab.q2))

# Collect quantile versions from the batches
source("../RPA/R/rpa.online.R")
qmat <- NULL
for (i in 1:length(batches)) {
  batch.file <- paste(batch.file.id, "-", i, ".RData", sep = "")
  load(batch.file)
  q <- get.probe.matrix(cels = batches[[i]], cdf = cdf, quantile.basis = quantile.basis, batch = batch)
  qmat <- cbind(qmat, q)
}

colnames(qmat) <- sapply(strsplit(colnames(qmat), "/"), function (x) {x[[8]]})

s <- colnames(qmat)[[1]]
cor(qmat[,s], qmat2[,s])


#system.time(qmat3 <- apply(pm(abatch.bg), 2, function(x) {quantile.basis[rank(x)]}))
#system.time(qmat4 <- set.quantiles(pm(abatch.bg), quantile.basis))

