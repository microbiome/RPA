set.seed(11122)
require(RPA)
fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}

################################################################

# Lukk Atlas
#cel.files <- list.celfiles("/homes/llahti/data/workdir/Lukk", full.names = T)
cel.files <- list.celfiles("CEL/", full.names = T)[1:12]
cdf <- "HGU133Plus2_Hs_ENSG" # cdf = NULL
batch.size <- 6

batch.file.id <- "batch.id"
hyperparameter.batches <- "batch.id.hyper"

bg.method <- "rma"
shuffle <- TRUE

###############################################################

message("Split CEL file list into batches")
batches <- get.batches(cel.files, batch.size, shuffle)

###############################################################

message("Calculating the basis for quantile normalization")

quantiles.done <- TRUE

if (!quantiles.done) {

  quantile.basis <- qnorm.basis.online(batches, bg.method = bg.method,
                                       cdf = cdf, save.batches =
                                       batch.file.id, batch.size =
                                       batch.size)
  
  save(quantile.basis, batches, file = "quantile.basis.RData")

} else {
    print("Load quantile basis and batch info")
    load("quantile.basis.RData")
}

###############################################################

hyperparameters.done <- TRUE

if (!hyperparameters.done) {
  
  hyps <- estimate.hyperparameters(batches = batches, cdf
                                   = cdf, quantile.basis =
                                   quantile.basis, bg.method =
                                   bg.method, load.batches =
                                   batch.file.id,
                                   save.hyperparameter.batches =
                                   hyperparameter.batches, mc.cores =
                                   1)

  save(hyps, file = "hyps.RData")

} else {
  print("Load hyperparameters")  
  load("hyps.RData")
}

###############################################################

summarization.done <- FALSE

if (!summarization.done) {

  # Final ExpressioSet object 
  message("Summarizing probesets")
  
  eset <- summarize.batches(variances = hyps$variances,
                            batches = batches, load.batches = batch.file.id,
                            mc.cores = 2, cdf = cdf,
                            quantile.basis = quantile.basis)
  
  # Arrange CEL files in the original order
  emat.rpa <- exprs(eset)[, cel.files]

  colnames(emat.rpa) <-
    sapply(strsplit(sapply(strsplit(colnames(emat.rpa), "/"),
    function(s){s[[length(s)]]}), "\\."), function(x){x[[1]]})
  
  save(emat.rpa, file = "Lukk-RPA.RData")

} else {
  load("Lukk-RPA.RData")
}

###################################################3

# Compare with RMA

abatch <- ReadAffy(filenames = cel.files)
if (!is.null(cdf)) {abatch@cdfName <- cdf}
emat.rma <- exprs(rma(abatch))
colnames(emat.rma) <-
    sapply(strsplit(sapply(strsplit(colnames(emat.rma), "/"),
    function(s){s[[length(s)]]}), "\\."), function(x){x[[1]]})
emat.rma <- emat.rma[, colnames(emat.rpa)]
if (!nrow(emat.rma) == nrow(emat.rpa)) {stop("Dimensions do not match")}
cors <- c()
for (i in 1:nrow(emat.rma)) {
    print(i) 
    cors[[i]] <- cor(emat.rpa[i, ], emat.rma[i,])
}

hist(cors)

######################################################
