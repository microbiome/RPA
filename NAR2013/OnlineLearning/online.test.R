require(RPA)
#fs <- list.files("../RPA/R", full.names = TRUE, pattern = ".R$")
#for (f in fs) {source(f)}

set.seed(11122) 

cel.files = list.celfiles("/share/mi/data/GSE3526/CEL/", full.names = T)[1:12]
                    sets = NULL
                  priors = list(alpha = 2, beta = 1)
                 epsilon = 1e-2 
                    cind = 1
                 verbose = FALSE
               bg.method = "rma"    
                 cdf <- "HGU133Plus2_Hs_ENSG" # cdf = NULL
              batch.size = 4
	  quantile.basis = NULL 
            save.batches = TRUE
  hyperparameter.batches = "batch.id.hyper"
                 shuffle = TRUE


  ###############################################################

  message("Split CEL file list into batches")
  batches <- get.batches(cel.files, batch.size, shuffle)

  ###############################################################

  message("Calculating the basis for quantile normalization")

  quantile.basis <- qnorm.basis.online(batches, bg.method, cdf, save.batches, batch.size)

  save(quantile.basis, file = "quantile.basis.RData")

  load("quantile.basis.RData")

  ###############################################################

  message("Estimating hyperparameters")
  hyps <- estimate.hyperparameters(sets, priors, batches, cdf, quantile.basis, bg.method, epsilon, cind, load.batches = batch.file.id, save.hyperparameter.batches = hyperparameter.batches, mc.cores = 2)

  save(batches, hyps, file = "hyps.RData")
  load("hyps.RData")

  ###############################################################

  # Final ExpressioSet object 
  message("Summarizing probesets")
  eset <- summarize.batches(sets, hyps$variances, batches, load.batches = batch.file.id, mc.cores = 2, cdf = cdf, quantile.basis = quantile.basis)

  save(eset, file = "eset.RData")
  load("eset.RData")

  ##################################################################

  # Arrange CEL files in the original order
  emat <- exprs(eset)
  emat <- emat[, cel.files]

  ##################################################################

  # RPA completed

  ##############################################

  # Compare to Original RPA and RMA

  # Remove path from CEL file names
  colnames(emat) <- unname(sapply(colnames(emat), function(s){unlist(strsplit(s, "/"))[[8]]}))

  # Compare
  abatch.tmp <- ReadAffy(filenames = cel.files)
  rpa <- RPA.pointestimate(abatch.tmp)

  eset.rma <- exprs(rma(abatch.tmp))
  eset.rpa <- exprs(rpa(abatch.tmp))
  s <- colnames(eset.rma)

  cors <- matrix(NA, nrow = nrow(eset.rma), ncol = 3)
  rownames(cors) <- rownames(eset.rma)
  colnames(cors) <- c("RMA.RPA", "RMA.RPAonline", "RPA.RPAonline")
  for (set in rownames(eset.rma)) {
    print(set)
    cors[set,] <- c(cor(eset.rma[set,s], eset.rpa[set,s]),
    cor(eset.rma[set,s], emat[set,s]),
    cor(eset.rpa[set,s], emat[set,s]))
  }

par(mfrow = c(3,1))
for (i in 1:ncol(cors)) {
  hist(cors[,i], 100, main = colnames(cors)[[i]])
}

  # Correlation between probe variance estimates from classical and online RPA
  X11(); plot(hyps$variances[[1]], rpa$sigma2[[1]]); abline(0,1)

k <- 1
print(cor(cbind(eset.rma[,s[[k]]], eset.rpa[,s[[k]]], emat[,s[[k]]])))





