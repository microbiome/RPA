# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

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
                  priors = list(alpha = 1, beta = 1),
                 epsilon = 1e-2,
                mc.cores = 1,
                 verbose = TRUE,                          
                 shuffle = TRUE,
              batch.size = 100, 
                 batches = NULL, 
          quantile.basis = NULL, 
            save.batches = TRUE,
	    save.batches.dir, 
	    keep.batch.files = FALSE, 
	    unique.run.identifier = NULL,
	    rseed = 23)

{

#fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names = T); for (f in fs) {source(f)}; cel.path = NULL; cel.files = cels; sets = NULL; cdf = NULL; bg.method = "rma"; priors = list(alpha = 1, beta = 1); epsilon = 1e-2; mc.cores = 4; verbose = TRUE; shuffle = TRUE; batch.size = bs; batches = NULL; quantile.basis = NULL; save.batches = TRUE; save.batches.dir = "."; keep.batch.files = FALSE; unique.run.identifier = NULL; rseed = 23

  ###############################################################

  # FIXME: make this working also with ready-made affybatches

  # Add a unique identifier for this RPA run
  if (is.null(unique.run.identifier)) {
    unique.run.identifier <- paste("RPA-run-id-", rnorm(1), sep = "")
  }

  if (!is.null(cel.path) && is.null(cel.files)) {
    message(paste("Preprocessing all CEL files from", cel.path))
    # Randomize the order to avoid biases in batch handling
    cel.files <- sample(list.celfiles(cel.path, full.names = T))
  }

  if (is.null(batch.size)) {
    if ( verbose ) { message("Determining batch size") }
    batch.size <- min(length(unlist(cel.files)), 100)
  }
  
  if (batch.size < 3) {
    warning("Minimum batch.size is 3. Setting batch.size = 3.")
    batch.size <- 3
  }

  if (is.null(batches)) {
    message("Split CEL file list into batches")
    batches <- get.batches(cel.files, batch.size, shuffle)
  }

  if (save.batches) {

    # Create the output directory if necessary
    if (length(dir(save.batches.dir)) == 0) { 
      message(paste("Creating directory", save.batches.dir))
      system(paste("mkdir ", save.batches.dir)) 
    }

    message(paste("Storing intermediate batches in directory", save.batches.dir, "with the identifier", unique.run.identifier))

    blf <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-batchlist.RData", sep = "")
    if ( verbose ) { message(paste("Saving batch list into file: ", blf)) }
    save(batches, file = blf)    
  }
  
  ###############################################################

  if (is.null(quantile.basis)) {
    message("Calculating the basis for quantile normalization")    
    quantile.basis <- qnorm.basis.online(batches, bg.method, cdf, save.batches = save.batches, batch.size, verbose = verbose, save.batches.dir = save.batches.dir, unique.run.identifier = unique.run.identifier)
  }
  if (save.batches) {
    quantile.file <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-quantiles.RData", sep = "")
    if (verbose) { message(paste("Saving quantile basis into file: ", quantile.file)) }
    save(quantile.basis, file = quantile.file)
  }

  ###############################################################
    
  hyper.parameters <- estimate.hyperparameters(sets, priors,
                                                  batches, cdf,
                                                  quantile.basis,
                                                  bg.method, epsilon,
                                                  load.batches = save.batches,
                                                  save.hyperparameter.batches = save.batches,
                                                  mc.cores = mc.cores,
                                                  verbose = verbose, 	    
						  save.batches.dir = save.batches.dir, 
						  unique.run.identifier = unique.run.identifier)
  
  if (save.batches) {
    hyper.file <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-hyperparameters.RData", sep = "")
    if ( verbose ) { message(paste("Saving hyperparameters into file: ", hyper.file)) }
    save(hyper.parameters, file = hyper.file)
  }

  ###############################################################

  if ( verbose ) { message("Estimating affinities..") }  
  affinities <- get.affinities(sets, hyper.parameters, cdf, mc.cores, quantile.basis, bg.method, batches, load.batches = TRUE, save.batches, save.batches.dir, unique.run.identifier, verbose)

  if ( verbose ) { message("Collect hyperparameters..") } # from batch files
  hyper.parameters <- get.probe.parameters(affinities, unique.run.identifier, save.batches.dir, mode = "list") 
  hyper.parameters.evolution <- collect.hyperparameters(batches, unique.run.identifier, save.batches.dir, batch.size)

  #################################################################

  # Given the affinities, calculate the final summary estimates
  message("Summarizing probesets")
  emat <- summarize.batches(sets = sets, 
       	  		    variances = hyper.parameters$tau2, 
			    batches = batches, 
			    load.batches = save.batches, 
		 	    mc.cores = mc.cores, 
			    cdf = cdf, 
			    bg.method = bg.method, 
			    quantile.basis = quantile.basis, 
			    verbose = verbose, 
			    save.batches.dir = save.batches.dir, 
			    unique.run.identifier = unique.run.identifier, 
			    save.batches = save.batches, 
			    affinities = affinities)

  ###############################################################

  if (!keep.batch.files) {
    message(paste("Removing the temporary batch files from directory", save.batches.dir, "with the identifier", unique.run.identifier))
    system(paste("rm ", save.batches.dir, "/", unique.run.identifier, "*", sep = ""))
  } else {
    message(paste("Keeping the temporary batch files in directory", save.batches.dir, "with the identifier", unique.run.identifier))
  }
  
  # Arrange CEL files in the original order and Coerce expression
  # values in the rpa object into an ExpressionSet object and return
  # expression set
  eset <- new("ExpressionSet", assayData = list(exprs = emat[, cel.files])) 

  ###############################################################

  # Store parameters
  params <- c(  cel.path = cel.path,
               cel.files = cel.files,
                    sets = sets,
                     cdf = cdf, 
               bg.method = bg.method,                              
                  priors = priors,
                 epsilon = epsilon,
                mc.cores = mc.cores,
                 verbose = verbose,                          
                 shuffle = shuffle,
              batch.size = batch.size, 
                 batches = batches, 
          quantile.basis = quantile.basis, 
            save.batches = save.batches,
	    save.batches.dir = save.batches.dir, 
	    keep.batch.files = keep.batch.files, 
	    unique.run.identifier = unique.run.identifier,
	    rseed = rseed)

  ###############################################################

  # Garbage collection
  gc()

  list(expressionSet = eset, hyper.parameters = hyper.parameters, hyper.parameters.evolution = hyper.parameters.evolution, params = params, sessionInfo = sessionInfo())

}



get.affinities <- function (sets, hyper.parameters, cdf, mc.cores, quantile.basis, bg.method, batches, load.batches = TRUE, save.batches, save.batches.dir, unique.run.identifier, verbose = TRUE) {

  if (verbose) { message("Estimating affinities") }

  # CEL files for this batch
  # (to speed up, estimate affinities based on the first batch only
  bi <- 1
  # Get background corrected, quantile normalized, and logged probe-level matrix
  batch <- NULL
  if (load.batches) {
    batch.file <- paste(save.batches.dir, "/", unique.run.identifier, "-", names(batches)[[bi]], ".RData", sep = "")
    if (verbose) { message(paste("Load preprocessed data for batch: ", batch.file)) }
    load(batch.file) # batch
  }

  # Get corresponding affinity estimates 
  # average affinity estimates over samples
  # this can be derived formally:
  # s_tj = a_r + d_t + mu_j + epsilon_tj; mu_j ~ N(0, tau2); ignore minimal variation from epsilon_tj; etc.
  # now, as a_t = a_r + d_t, we have a_t ~ s_tj - mu_j ~ N(s_tj, tau2_j) -> ML estimate of a_t is given by 
  # weighted average of s_tj, weighted by tau2_j's
  batch.cels <- batches[[bi]]
  set.inds <- get.set.inds(batch.cels[1:2], cdf, sets)
  if ( is.null(sets) ) { sets <- names(set.inds) }
    
  # Get probes x samples matrices for each probeset
  # No need to remove the reference sample for d.update here in the summarization step!
  # Since probe-specific variance is now known (from the estimation above), the probeset-level signal
  # estimate is obtained as a weighted sum of the probes, weighted by the probe-specific variances
  if ( verbose ) { message("Extract probe-level data") }       
  q <- get.probe.matrix(cels = batch.cels, cdf = cdf, quantile.basis = quantile.basis, bg.method = bg.method, batch = batch, verbose = verbose)
  q <- mclapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) }, mc.cores = mc.cores) 
  names(q) <- sets

  # Get the summary estimate for probeset using the posterior variances
  # (>25-fold speedup with sapply)
  variances <- hyper.parameters$tau2
  em <- t(sapply(mclapply(sets, function(set) {d.update.fast(q[[set]], variances[[set]])}, mc.cores = mc.cores), identity))
  rownames(em) <- sets
  colnames(em) <- batch.cels

  # Given probe-level observations and summary estimates, calculate the corresponding affinities
  affinities <- lapply(sets, function (set) {estimate.affinities(q[[set]], em[set, batch.cels])})
  names(affinities) <- sets 

  affinities

}

summarize.batches <- function (sets = NULL, variances, batches, load.batches = FALSE, mc.cores = 1, cdf = NULL, bg.method = "rma", normalization.method = "quantiles", verbose = TRUE, quantile.basis, save.batches.dir = ".", unique.run.identifier = NULL, save.batches = FALSE, affinities) {

  #sets = sets; variances = hyper.parameters$variances; batches = batches; load.batches = save.batches; mc.cores = mc.cores; cdf = cdf; bg.method = bg.method; quantile.basis = quantile.basis; verbose = verbose; normalization.method = "quantiles"
     
  # FIXME: remove normalization method from here as unnecessary?
  if ( verbose ) { message("Pick PM indices") }
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
    if (load.batches) {
      batch.file <- paste(save.batches.dir, "/", unique.run.identifier, "-", names(batches)[[i]], ".RData", sep = "")
      if (verbose) { message(paste("Load preprocessed data for this batch from: ", batch.file)) }
      load(batch.file) # batch
    }
    
    # Get probes x samples matrices for each probeset
    # No need to remove the reference sample for d.update here in the summarization step!
    # Since probe-specific variance is now known (from the estimation above), the probeset-level signal
    # estimate is obtained as a weighted sum of the probes, weighted by the probe-specific variances
    if ( verbose ) { message("Extract probe-level data") }       
    q <- get.probe.matrix(cels = batch.cels, cdf = cdf, quantile.basis = quantile.basis, bg.method = bg.method, batch = batch, verbose = verbose)
    q <- mclapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) }, mc.cores = mc.cores) 
    names(q) <- sets

    # Get summary estimates for each probeset with posterior variances
    # Remove the estimated affinity effects from the raw signals at each probe before summarization
    # (>25-fold speedup obtained with sapply)
    if ( verbose ) { message("Summarizing..") }       
    emat[sets, batch.cels] <- t(sapply(mclapply(sets, function(set) {d.update.fast(q[[set]] - affinities[[set]], variances[[set]])}, mc.cores = mc.cores), identity))

  }

  emat
}  







get.probe.matrix <- function (cels, cdf = NULL, quantile.basis, bg.method = "rma", normalization.method = "quantiles", batch = NULL, verbose = TRUE) {
  
  if (!is.null(batch)) {

      # Assuming that the bg correction + quantile normalization have
      # been already calculated for quantile.basis, which is here
      # simply allocated for each array
    
      if (verbose) { message("Set quantile data on each array") }
      q <- apply(batch, 2, function (o) {quantile.basis[o]}) 
      if (verbose) { message("...Done.") }
      
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
        if (verbose) { message("Normalizing") }
        pm(abatch) <- set.quantiles(pm(abatch), quantile.basis)
      } 
  
      if (verbose) {message("Taking logarithm")}
      q <- log2(pm(abatch))

  }

  q

}

