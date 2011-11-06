
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
	       normalization.method = "quantiles",
                     cdf = NULL, 
              batch.size = 10, 
	  quantile.basis = NULL) 
{

  warning("rpa.online is an experimental version")

  if (batch.size < 3) {
    warning("Minimum batch.size is 3. Setting batch.size = 3.")
    batch.size <- 3
  }

  if (is.null(cel.files) && !is.null(cel.path)) {
    cel.files <- list.celfiles(cel.path, full.names = TRUE)
  }

  #################################################################
  
  # NOTE: list CEL files in random order to avoid biases!
  # FIXME: add this
  cel.files.sampled <- cel.files # sample(cel.files)

  # Get probe position indices
  abatch <- ReadAffy(filenames = cel.files[1:2], compress=getOption("BioC")$affy$compress.cel) 
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
  
  ###############################################################

  if (normalization.method == "quantiles") {
    # Online-estimation of the basis for quantile normalization
    message("Calculating the basis for quantile normalization")
    quantile.basis <- quantile.basis.online(cel.files.sampled, bg.method, batch.size, cdf)
  }
  
  ###############################################################
  
  # Split CEL file list into batches
  # TODO: list CEL files in random order to avoid biases!
  # TODO: Define batching already before quantile.basis calculation?
  batches <- get.batches(cel.files.sampled, batch.size)
 
  ###############################################################

  hyps <- estimate.hyperparameters(priors, sets, set.inds, batches, cdf, quantile.basis, bg.method, normalization.method, epsilon)
  alpha <- hyps$alpha  
  betas <- hyps$betas

  ###############################################################

  # Get final estimated variances for each probeset based on hyperparameter posteriors
  variances <- lapply(betas, function (beta) {beta/alpha})
  names(variances) <- names(betas) 
  
  # Now hyperparameters for probe variances have been estimated
  # Get probeset-level signal estimate by online sweep with the fixed variances
  emat <- array(NA, dim = c(length(sets), length(cel.files.sampled)))
  rownames(emat) <- sets
  colnames(emat) <- cel.files.sampled # sapply(strsplit(cel.files, "/"), function (x) {x[[length(x)]]})

  for (i in 1:length(batches)) {

    message(paste("Summarizing batch", i, "/", length(batches)))

    # CEL file IDs for this batch
    batch.cels <- batches[[i]]
    
    # Get background corrected, quantile normalized, and logged probe-level matrix
    # Do NOT calculate probe-level diff.exp here any more, only needed in variance estimation.
    q <- get.probe.matrix(cels = batch.cels, cdf, quantile.basis, bg.method, normalization.method)

    # Get probes x samples matrices for each probeset
    q <- lapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) })
    names(q) <- sets

    for (set in sets) {
      # print(which(set == sets)/length(sets))    
      # Get summary estimate using the posterior variance
      emat[set, batch.cels] <- d.update.fast(q[[set]], variances[[set]])
    }
  }

  # Return the arrays into original order
  emat <- emat[, cel.files]
  
  # Remove path from CEL names for compactness
  # colnames(emat) <- sapply(strsplit(colnames(emat), "/"), function (x) { x[[length(x)]] })

  # Coerce expression values in the rpa object into an ExpressionSet object
  # and return expression set
  new("ExpressionSet", assayData = list(exprs = emat)) 
  
}

summarize.batches <- function (sets, variances, batches, load.batches, mc.cores = 1) {

  cel.files <- unlist(batches)		  

  # Now hyperparameters for probe variances have been estimated
  # Get probeset-level signal estimate by online sweep with the fixed variances
  emat <- array(NA, dim = c(length(sets), length(cel.files)))
  rownames(emat) <- sets
  colnames(emat) <- cel.files # sapply(strsplit(cel.files, "/"), function (x) {x[[length(x)]]})

  for (i in 1:length(batches)) {

    message(paste("Summarizing batch", i, "/", length(batches)))

    # CEL files for this batch
    batch.cels <- batches[[i]] 

    # Get background corrected, quantile normalized, and logged probe-level matrix

    if (!is.null(load.batches)) {
      batch.file <- paste(load.batches, "-", i, ".RData", sep = "")
    }
      
    q <- get.probe.matrix(cels = batch.cels, cdf, quantile.basis, bg.method, normalization.method, batch.file, verbose = verbose)

    # Get probes x samples matrices for each probeset
    # No need to remove the reference sample for d.update here
    # FIXME with mc.cores = 1 mclapply -> lapply, no need to separate
    q <- mclapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) }, mc.cores = mc.cores) 
    names(q) <- sets

    for (set in sets) {
      # Get summary estimate for each probeset using the posterior variance
      emat[set, batch.cels] <- d.update.fast(q[[set]], variances[[set]])
    }
  }

  # Remove path from CEL names for compactness
  # colnames(emat) <- sapply(strsplit(colnames(emat), "/"), function (x) { x[[length(x)]] })

  # Coerce expression values in the rpa object into an ExpressionSet object
  # and return expression set
  new("ExpressionSet", assayData = list(exprs = emat)) 
}  

get.probe.matrix <- function (cels, cdf, quantile.basis, bg.method = "rma", normalization.method = "quantiles", batch.file = NULL, verbose = TRUE) {

  if (!is.null(batch.file)) {
      if (verbose) {message("Load precalculated batch")}
      # (gives probe ordering for quantile normalization)
      # Intermediate saves should remarkably save time here in calculations
      #load(paste(load.batches, "-", i, ".RData", sep = "")) # batch: apply(pm(abatch), 2, order)
      load(batch.file)
      
      if (verbose) {message("Set quantile data on each array")}
      q <- apply(batch, 2, function (o) {quantile.basis[o]}) 
      #q <- set.quantiles(pm(abatch.bg), qb)
    } else {
  
  # Get background corrected, quantile normalized, and logged probe-level matrix
  # The bg + quant are from quantile.basis  

  # Getting affybatch
  abatch <- ReadAffy(filenames = cels, compress=getOption("BioC")$affy$compress.cel) 
  
  # Set alternative CDF environment if given
  if (!is.null(cdf)) { 
    message("Setting altCDF..")
    abatch@cdfName <- cdf   
  } else {
    message("Default CDF being used..")
  }
  
  if (bg.method == "none") {
    message("Background correction skipped..")
  } else {
    message("Background correcting...")
    abatch <- bg.correct(abatch, bg.method, destructive = TRUE)
  }

  # Normalize by forcing the pre-calculated quantile basis
  if (normalization.method == "none") {
    message("Normalization skipped..")
  } else if (normalization.method == "quantiles") {
    message("Normalizing...")
    warning("Remove this saving from release")
    save(abatch, quantile.basis, file = "tmp.RData")
    pm(abatch) <- set.quantiles(pm(abatch), quantile.basis)
  } 
  
  # Log transformation
  q <- log2(pm(abatch))

}

  q

}

