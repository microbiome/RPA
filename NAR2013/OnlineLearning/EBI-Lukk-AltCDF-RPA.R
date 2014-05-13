set.seed(11122)
require(RPA)

#fs <- list.files("../RPA/R", full.names = TRUE, pattern = ".R$")
#for (f in fs) {source(f)}

################################################################

quantiles.done <- FALSE
hyperparameters.done <- FALSE
summarization.done <- FALSE

# Lukk Atlas                                        
#cel.files <- list.celfiles("/homes/llahti/data/workdir/Lukk", full.names = T)
#cel.files <- list.celfiles("/share/mi/data/GSE3526/CEL/", full.names = T)[1:50]
cel.files <- list.celfiles("CEL/", full.names = T)

#cdf <- "HGU133Plus2_Hs_ENSG" # cdf = NULL

batch.size <- 200

save.batches <- TRUE
bg.method <- "rma"
shuffle <- TRUE

###############################################################

message("Split CEL file list into batches")
batches <- get.batches(cel.files, batch.size, shuffle)

###############################################################

message("Calculating the basis for quantile normalization")

if (!quantiles.done) {

  quantile.basis <- qnorm.basis.online(batches, bg.method = bg.method,
                                       cdf = cdf, save.batches =
                                       save.batches, batch.size =
                                       batch.size)

  save(quantile.basis, batches, file = "quantile.basis.RData")

} else {
    print("Load quantile basis and batch info")
    load("quantile.basis.RData")
}

###############################################################

if (!hyperparameters.done) {
  
  hyps <- estimate.hyperparameters(batches = batches, cdf = cdf,
                                   quantile.basis = quantile.basis,
                                   bg.method = bg.method, load.batches
                                   = save.batches,
                                   save.hyperparameter.batches =
                                   save.batches, mc.cores = 1)

  save(hyps, file = "hyps.RData")

} else {
  print("Load hyperparameters")  
  load("hyps.RData")
}

###############################################################


if (!summarization.done) {

  # Final ExpressioSet object 
  message("Summarizing probesets")
  
  emat <- summarize.batches(variances = hyps$variances,
                            batches = batches, load.batches = save.batches,
                            mc.cores = 2, cdf = cdf,
                            quantile.basis = quantile.basis)
  
  # Arrange CEL files in the original order
  emat <- emat[, cel.files]

  colnames(emat) <-
    sapply(strsplit(sapply(strsplit(colnames(emat), "/"),
                           function(s){s[[length(s)]]}), "\\."), function(x){x[[1]]})
  
  save(emat, file = "Lukk-RPA.RData")

} else {
  load("Lukk-RPA.RData")
}

###################################################3

# Load Lukk RMA

check <- TRUE
if (check) {

  # Compare to Original RPA and RMA

  # Compare
  abatch.tmp <- ReadAffy(filenames = cel.files)
  if (!is.null(cdf)) { abatch.tmp@cdfName <- cdf }
  
  # rpa <- RPA.pointestimate(abatch.tmp, cdf = cdf)

  eset.rma <- exprs(rma(abatch.tmp))
  eset.rpa <- exprs(rpa(abatch.tmp, cdf = cdf))

  colnames(eset.rma) <-
    sapply(strsplit(sapply(strsplit(colnames(eset.rma), "/"),
                           function(s){s[[length(s)]]}), "\\."), function(x){x[[1]]})

  colnames(eset.rpa) <-
    sapply(strsplit(sapply(strsplit(colnames(eset.rpa), "/"),
                           function(s){s[[length(s)]]}), "\\."), function(x){x[[1]]})  


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

######################################################
#source("lukkrma.R")
