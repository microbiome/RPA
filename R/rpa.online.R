
# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  


summarize.batches <- function (sets, variances, batches, load.batches = NULL, save.batches = NULL, mc.cores = 1, cdf = NULL, bg.method = NULL, normalization.method = NULL, verbose = FALSE, quantile.basis) {

#load.batches = batch.file.id; save.batches = NULL; mc.cores = 2; cdf = NULL; bg.method = NULL; normalization.method = NULL; verbose = FALSE

  
  # FIXME: remove normalization method from here as unnecessary?
  cel.files <- unlist(batches)		  

  # Now hyperparameters for probe variances have been estimated
  # Get probeset-level signal estimate by online sweep with the fixed variances
  emat <- array(NA, dim = c(length(sets), length(cel.files)))
  rownames(emat) <- sets
  colnames(emat) <- cel.files # sapply(strsplit(cel.files, "/"), function (x) {x[[length(x)]]})

  if (verbose) {message("Pick PM indices")}
  set.inds <- get.set.inds(batches[[1]][1:2], cdf)
  
  for (i in 1:length(batches)) {

    message(paste("Summarizing batch", i, "/", length(batches)))

    # CEL files for this batch
    batch.cels <- batches[[i]] 

    # Get background corrected, quantile normalized, and logged probe-level matrix

    if (!is.null(load.batches)) {
      batch.file <- paste(load.batches, "-", i, ".RData", sep = "")
      load(batch.file) # batch
    } else {
      batch <- NULL
    }
      
    q <- get.probe.matrix(cels = batch.cels, cdf, quantile.basis, bg.method, normalization.method, batch, verbose = verbose)

    # Get probes x samples matrices for each probeset
    # No need to remove the reference sample for d.update here
    # FIXME with mc.cores = 1 mclapply -> lapply, no need to separate
    
    q <- mclapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) }, mc.cores = mc.cores) 
    names(q) <- sets

    # FIXME: use mclapply
    for (set in sets) {
      # Get summary estimate for each probeset using the posterior variance
      emat[set, batch.cels] <- d.update.fast(q[[set]], variances[[set]])
    }

    # emat[, batch.cels] <- as.matrix(mclapply(set, function(set) {d.update.fast(q[[set]], variances[[set]])}))
    #
    #

    
  }

  # Remove path from CEL names for compactness
  # colnames(emat) <- sapply(strsplit(colnames(emat), "/"), function (x) { x[[length(x)]] })

  # Coerce expression values in the rpa object into an ExpressionSet object
  # and return expression set
  new("ExpressionSet", assayData = list(exprs = emat)) 
}  

get.probe.matrix <- function (cels, cdf = NULL, quantile.basis, bg.method = "rma", normalization.method = "quantiles", batch = NULL, verbose = TRUE) {

#  cels = batch.cels

  
  if (!is.null(batch)) {
      
      if (verbose) {message("Set quantile data on each array")}
      q <- apply(batch, 2, function (o) {quantile.basis[o]}) 
      if (verbose) {message("...Done.")}
      
    } else {
  
      # Get background corrected, quantile normalized, and logged
      # probe-level matrix The bg + quant are from quantile.basis

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
        pm(abatch) <- set.quantiles(pm(abatch), quantile.basis)
      } 
  
      if (verbose) {message("Taking logarithm")}
      q <- log2(pm(abatch))

    }

  q

}

