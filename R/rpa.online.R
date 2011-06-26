# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  
rpa.online <- function (
 	        cel.path = NULL,
               cel.files = NULL, 
                    sets = NULL,
                  priors = list(alpha = 2, beta = 1),
                 epsilon = 1e-2, 
                    cind = 1,
                 verbose = FALSE,
               bg.method = "rma",
                     cdf = NULL, 
              batch.size = 10) 
{

  warning("rpa.online is an experimental version")

  if (batch.size < 2) {
     warning("Minimum batch.size is 2. Setting batch.size = 2.")
     batch.size <- 2
  }

  if (is.null(cel.files) && !is.null(cel.path)) {
       cel.files <- list.celfiles(cel.path, full.names = TRUE)
  }

  # NOTE: list CEL files in random order to avoid biases!
  cel.files <- sample(cel.files)

  ###############################################################

  # Online-estimation of the basis for quantile normalization
  message("Calculating the basis for quantile normalization")
  quantile.basis <- quantile.basis.online(cel.files, bg.method, batch.size, cdf)

  ###############################################################
  
  # Split CEL file list into batches
  # NOTE: list CEL files in random order to avoid biases!
  batches <- get.batches(cel.files, batch.size)

  # Get probe position indices
  # Getting first affybatch
  abatch <- ReadAffy(filenames = batches[[1]], compress=getOption("BioC")$affy$compress.cel) 
  # Set alternative CDF environment if given
  if (!is.null(cdf)) { abatch@cdfName <- cdf }    

  # Check names and number for the investigated probesets
  # if my.sets not specified, take all sets in abatch
  if (is.null(sets)) { sets <- geneNames(abatch) } 
  Nsets <- length( sets ) 

  # Retrieve probe positions
  # The indices correspond to the output of pm()
  pN <- probeNames(abatch)
  set.inds <- split(1:length(pN), pN)[sets] # pNList

  #################################################################

  # Hyperparameter estimation through batches

  # Initialize hyperparameters
  # Note: alpha is scalar and same for all probesets alpha <- alpha + N/2 at each batch
  alpha <- priors$alpha # initialize 
  betas <- lapply(sets, function (set) { 
    rep(priors$beta, length(set.inds[[set]]))
  })
  names(betas) <- sets    

  for (i in 1:length(batches)) {

    message(paste("Updating hyperparameters on batch", i, "/", length(batches)))

    # Get background corrected, quantile normalized, and logged probe-level matrix
    q <- get.probe.matrix(cels = batches[[i]], cdf, quantile.basis)

    # Get probes x samples matrix of probe-wise fold-changes
    q <- matrix(q[, -cind] - q[, cind], nrow(q))

    T <- ncol(q) # Number of arrays expect reference

    # Get probes x samples matrices for each probeset
    q <- lapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) })
    names(q) <- sets

    # Update variance for each probeset
    alpha.old <- alpha 
    s2s <- lapply(sets, function (set) {
      s2.update(q[[set]], alpha.old, betas[[set]], s2.init = betas[[set]]/alpha.old, th = 1e-2)
    }) # FIXME move conv. param. to arguments
    names(s2s) <- sets

    # Update alpha, beta (variance = beta/alpha at mode with large T)
    alpha <- update.alpha(T, alpha)
    betas <- lapply(s2s, function (s2) {s2 * alpha})

  }

  # Now hyperparameters for probe variances have been estimated
  # Get probeset-level signal estimate by online sweep with the fixed variances

  emat <- array(NA, dim = c(length(sets), length(cel.files)))
  rownames(emat) <- sets
  colnames(emat) <- sapply(strsplit(cel.files, "/"), function (x) {x[[length(x)]]})

  for (i in 1:length(batches)) {

    message(paste("Summarizing batch", i, "/", length(batches)))

    batch.cels <- batches[[i]]

    # Get background corrected, quantile normalized, and logged probe-level matrix
    # Do NOT calculate probe-level diff.exp here any more, only needed in variance estimation.
    q <- get.probe.matrix(cels = batch.cels, cdf, quantile.basis)

    # Get probes x samples matrices for each probeset
    q <- lapply(set.inds, function (pmis) {matrix(q[pmis,], length(pmis))})
    names(q) <- sets

    # Get estimated variances for each probeset based on hyperparameter posteriors
    variances <- lapply(betas, function (beta) {beta/alpha})
    names(variances) <- names(betas) 

    for (set in sets) {
      # print(which(set == sets)/length(sets))    
      # Get summary estimate using the posterior variance
      cel.names <- sapply(strsplit(batch.cels, "/"), function (x) { x[[length(x)]] })
      emat[set, cel.names] <- d.update.fast.c(q[[set]], variances[[set]])
      # d.update.fast(q[[set]], variances[[set]])
    }
  }

  # Coerce expression values in the rpa object into an ExpressionSet object
  # and return expression set
  new("ExpressionSet", assayData = list(exprs = emat)) 
  
}


get.probe.matrix <- function (cels, cdf, quantile.basis) {

  # Get background corrected, quantile normalized, and logged probe-level matrix
  # The bg + quant are from quantile.basis  

  # Getting affybatch
  abatch <- ReadAffy(filenames = cels, compress=getOption("BioC")$affy$compress.cel) 

  # Set alternative CDF environment if given
  if (!is.null(cdf)) { abatch@cdfName <- cdf }    

  # Normalize by forcing the pre-calculated quantile basis
  pm(abatch) <- set.quantiles(pm(abatch), quantile.basis)

  # Log transformation
  log2(pm(abatch))

}