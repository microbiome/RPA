
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
                     cdf = NULL, 
               bg.method = "rma",                              
                  priors = list(alpha = 2, beta = 1),
                 epsilon = 1e-2,

                    cind = 1,
                mc.cores = 1,
                 verbose = TRUE,                          
                 shuffle = TRUE,
                          
              batch.size = 10, 
                 batches = NULL, # user-defined CEL file batches                          
           batch.file.id = NULL, # will be augmented with batch number and .RData suffix
      batch.list.file.id = NULL,                          
                          
          quantile.basis = NULL, # pre-calculated basis for quantile normalization
        quantile.file.id = NULL, # will be augmented with .RData suffix
                          
        hyper.parameters = NULL, 
 hyperparameter.batch.id = NULL # will be augmented with batch number and .RData suffix
                          )     # except on the last round when all batches have been considered
                                # then no batch number used
{

  ###############################################################

  if (is.null(batches)) {
    message("Split CEL file list into batches")
    batches <- get.batches(cel.files, batch.size, shuffle)
  }
  if (!is.null(batch.list.file.id)) {
    if (verbose) {message(paste("Saving batch list into file: ", batch.list.file.id))}
    save(batches, file = batch.list.file.id)    
  }
  
  
  ###############################################################

  if (is.null(quantile.basis)) {
    message("Calculating the basis for quantile normalization")    
    quantile.basis <- qnorm.basis.online(batches, bg.method, cdf, batch.file.id, batch.size, verbose = verbose)
  }
  if (!is.null(quantile.file.id)) {
    quantile.file <- paste(quantile.file.id, ".RData", sep = "")
    if (verbose) {message(paste("Saving quantile basis into file: ", quantile.file))}
    save(quantile.basis, file = quantile.file)
  }

  ###############################################################

  if (is.null(hyper.parameters)) {
    message("Estimating hyperparameters")
    hyper.parameters <- estimate.hyperparameters(sets, priors, batches, cdf, quantile.basis, 
                                                  bg.method, epsilon, cind, load.batches = batch.file.id, 
						  save.hyperparameter.batches = hyperparameter.batch.id, 
						  mc.cores = mc.cores, verbose = verbose)
  }
  
  if (!is.null(hyperparameter.batch.id)) {
    hyper.file <- paste(hyperparameter.batch.id, ".RData", sep = "")
    if (verbose) {message(paste("Saving hyperparameters into file: ", hyper.file))}
    save(hyper.parameters, file = hyper.file)
  }

  ###############################################################

  # Final ExpressioSet object 
  message("Summarizing probesets")
  eset <- summarize.batches(sets, hyper.parameters$variances, batches, load.batches = batch.file.id, mc.cores = mc.cores, cdf = cdf, quantile.basis = quantile.basis, verbose = verbose)
  # save(eset, file = "eset.RData")
  
  ##################################################################

  # Arrange CEL files in the original order and Coerce expression
  # values in the rpa object into an ExpressionSet object and return
  # expression set
  new("ExpressionSet", assayData = list(exprs = exprs(eset)[, cel.files])) 

}


summarize.batches <- function (sets = NULL, variances, batches, load.batches = NULL, save.batches = NULL, mc.cores = 1, cdf = NULL, bg.method = NULL, normalization.method = NULL, verbose = FALSE, quantile.basis) {

  # FIXME: remove normalization method from here as unnecessary?
  if (verbose) {message("Pick PM indices")}
  set.inds <- get.set.inds(batches[[1]][1:2], cdf, sets)
  if (is.null(sets)) {sets <- names(set.inds)}
  
  # Now hyperparameters for probe variances have been estimated
  # Get probeset-level signal estimate by online sweep with the fixed variances
  cel.files <- unlist(batches)
  emat <- array(NA, dim = c(length(sets), length(cel.files)))
  rownames(emat) <- sets
  colnames(emat) <- cel.files # sapply(strsplit(cel.files, "/"), function (x) {x[[length(x)]]})
  
  for (i in 1:length(batches)) {

    message(paste("Summarizing batch", i, "/", length(batches)))

    # CEL files for this batch
    batch.cels <- batches[[i]] 

    # Get background corrected, quantile normalized, and logged probe-level matrix
    batch <- NULL
    if (!is.null(load.batches)) {
      batch.file <- paste(load.batches, "-", i, ".RData", sep = "")
      if (verbose) {message(paste("Load preprocessed data for this batch from: ", batch.file))}
      load(batch.file) # batch
    }

    # Get probes x samples matrices for each probeset
    # No need to remove the reference sample for d.update in the summarization step!
    if (verbose) {message("Extract probe-level data")}      
    q <- get.probe.matrix(batch.cels, cdf, quantile.basis, bg.method, normalization.method, batch, verbose = verbose)    
    q <- mclapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) }, mc.cores = mc.cores) 
    names(q) <- sets

    # Get summary estimate for each probeset using the posterior variance
    emat[sets, batch.cels] <- t(sapply(mclapply(sets, function(set) {d.update.fast(q[[set]], variances[[set]])}, mc.cores = mc.cores), identity))
    # Check 11.11.2011 - more than 25-fold speedup
    # system.time(for (set in sets) {emat[set, batch.cels] <- d.update.fast(q[[set]], variances[[set]])})    

  }

  # Coerce expression values in the rpa object into an ExpressionSet object
  # and return expression set
  new("ExpressionSet", assayData = list(exprs = emat)) 
}  

get.probe.matrix <- function (cels, cdf = NULL, quantile.basis, bg.method = "rma", normalization.method = "quantiles", batch = NULL, verbose = TRUE) {
  
  if (!is.null(batch)) {

      # Assuming that the bg correction + quantile normalization have
      # been already calculated for quantile.basis, which is here
      # simply allocated for each array
    
      if (verbose) {message("Set quantile data on each array")}
      q <- apply(batch, 2, function (o) {quantile.basis[o]}) 
      if (verbose) {message("...Done.")}
      
  } else {
  
      # Calculate background corrected, quantile normalized, and
      # logged probe-level matrix from the CEL files if precalculated
      # batches are not available.

      # Getting affybatch
      abatch <- ReadAffy(filenames = cels, compress=getOption("BioC")$affy$compress.cel) 
  
      # Set alternative CDF environment if given
      if (!is.null(cdf)) { 
        message("Setting altCDF..")
        abatch@cdfName <- cdf   
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

