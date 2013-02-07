# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' rpa.online
#' RPA-online for preprocessing very large expression data sets.
#'
#' @param cel.path Path to CEL file directory
#' @param cel.files List of CEL files to preprocess
#' @param sets Probesets for which RPA will be computed
#' @param cdf Specify an alternative CDF environment
#' @param bg.method Specify background correction method. See bgcorrect.methods() for options.
#' @param priors  An 'rpa.priors' object. Can be used to set user-specified priors for the model parameters. Not used tau2.method = "var". The prior parameters alpha and beta are prior parameters for inverse Gamma distribution of probe-specific variances. Noninformative prior is obtained with alpha, beta -> 0.  Not used with tau2.method 'var'. Scalar alpha and beta specify an identical inverse Gamma prior for all probes, which regularizes the solution. Can be also specified as lists, each element corresponding to one probeset.
#' @param epsilon  Convergence tolerance. The iteration is deemed converged when the change in all parameters is < epsilon.
#' @param mc.cores Number of cores for parallel computation
#' @param verbose Print progress information during computation
#' @param shuffle Form random batches  
#' @param batch.size Batch size for online mode (rpa.online); the complete list of CEL files will be preprocessed in batches with this size using Bayesian online-updates for probe-specific parameters.
#' @param batches User-defined CEL file batches
#' @param quantile.basis Pre-calculated basis for quantile normalization
#' @param save.batches.dir Output directory for temporary batch saves.
#' @param keep.batch.files Logical. Keep (TRUE) or remove (FALSE) the batch files after preprocessing.
#' @param unique.run.identifier Define identifier for this run for naming the temporary batch files. By default, a random id is generated.
#' @param rseed Random seed.
#'
#' @details rpa.online is used to preprocess very large expression data collections based on a Bayesian hyperparameter update procedure. Returns an expressionSet object preprocessed with RPA. Gives an estimate of the probeset-level mean parameter d of the RPA model, and returns these in an expressionSet object. The CEL files are handled in batches to obtain Bayesian updates for probe-specific hyperpriors; after sweeping through the database in batches the results are combined. The online mode is useful for preprocessing very large expression data sets where ordinary preprocessing algorithms fail, without compromises in modelling stage.
#'
#' @return An instance of the 'expressionSet' class.
#'
#' @seealso rpa, AffyBatch, ExpressionSet 
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # eset <- rpa.online(cel.file.path) 
#' @keywords methods

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
	    save.batches.dir = ".", 
	    keep.batch.files = FALSE, 
	    unique.run.identifier = paste("RPA-run-id-", rnorm(1), sep = ""),
	    rseed = 23)

{

  # cel.path = "CEL"; cel.files = NULL; sets = NULL; cdf = NULL;  bg.method = "rma"; priors = list(alpha = 1, beta = 1); epsilon = 1e-2; mc.cores = 1; verbose = TRUE; shuffle = TRUE; batches = NULL; quantile.basis = NULL; save.batches.dir = "."; keep.batch.files = FALSE; unique.run.identifier = paste("RPA-run-id-", rnorm(1), sep = ""); rseed = 23; batch.size = Inf; save.batches = FALSE

  # Always store the batches temporarily
  save.batches <- TRUE

  # Get CEL file list
  if (!is.null(cel.path) && is.null(cel.files)) {
    message(paste("Preprocessing all CEL files from", cel.path))
    cel.files <- list.celfiles(cel.path, full.names = T)
  }

  # Define batches
  if ( is.null(batches) ) {
    message("Split CEL file list into batches")
    batches <- get.batches(cel.files, batch.size, shuffle)
  }

  # Optionally save batch list
  blf <- saving.batchlist(batches, save.batches, save.batches.dir, unique.run.identifier, verbose)
  
  # ---------------------------------------------------------------------------

  # CALCULATE THE BASIS FOR QUANTILE NORMALIZATION
  if (is.null(quantile.basis)) {
    quantile.basis <- qnorm.basis.online(batches, bg.method, cdf, save.batches = save.batches, verbose = verbose, save.batches.dir = save.batches.dir, unique.run.identifier = unique.run.identifier)
  }

  # Optionally save quantile basis
  qf <- saving.quantile.basis(quantile.basis, save.batches, save.batches.dir, unique.run.identifier, verbose)

  # ---------------------------------------------------------------------------------

  # HYPERPARAMETER ESTIMATION

  if (verbose) { message("Get probeset indices") }
  set.inds <- get.set.inds(batches[[1]][1:2], cdf, sets)

  hyper.parameters <- estimate.hyperparameters(sets, priors,
                                                  batches, cdf,
                                                  quantile.basis,
                                                  bg.method, epsilon,
                                                  load.batches = save.batches,
                                                  save.hyperparameter.batches = save.batches,
                                                  mc.cores = mc.cores,
                                                  verbose = verbose, 	    
						  save.batches.dir = save.batches.dir, 
						  unique.run.identifier = unique.run.identifier, 
						  set.inds = set.inds)

  hf <- saving.hyperparameters(hyper.parameters, save.batches, save.batches.dir, unique.run.identifier, verbose)

  # -----------------------------------------------------------------------------------

  # AFFINITY ESTIMATION
  affinities <- get.affinities(sets, hyper.parameters, cdf, mc.cores, quantile.basis, bg.method, batches, load.batches = save.batches, save.batches, save.batches.dir, unique.run.identifier, verbose, set.inds)

  # Add affinities into hyperparameter list
  hyper.parameters$mu <- affinities

  # ----------------------------------------------------------------------------------------------
  
  # COLLECTING HYPERPARAMETERS ACROSS ALL BATCHES TO INVESTIGATE HYPERPARAMETER EVOLUTION
  hyper.parameters.evolution <- collect.hyperparameters(batches, unique.run.identifier, save.batches.dir, save.batches)

  # -----------------------------------------------------------------------------------

  # PROBESET SUMMARIZATION

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
			    affinities = affinities, set.inds = set.inds)

  # Arrange CEL files in the original order and Coerce expression
  # values in the rpa object into an ExpressionSet object and return
  # expression set object
  eset <- new("ExpressionSet", assayData = list(exprs = emat[, cel.files])) 

  # ---------------------------------------------------------------------------

  # Remove temporary batch files and run garbage collection
  tmp <- batch.cleanup(keep.batch.files, unique.run.identifier, save.batches.dir)

  # -------------------------------------------------------------------------------------

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
	      batch.size = length(batches[[1]]),
                 batches = batches, 
          quantile.basis = quantile.basis, 
            save.batches = save.batches,
        save.batches.dir = save.batches.dir, 
	keep.batch.files = keep.batch.files, 
   unique.run.identifier = unique.run.identifier,
	           rseed = rseed, 
	     sessionInfo = sessionInfo())

  list(expressionSet = eset, hyper.parameters = hyper.parameters, hyper.parameters.evolution = hyper.parameters.evolution, params = params)

}



batch.cleanup <- function (keep.batch.files, unique.run.identifier, save.batches.dir) {

  if (!keep.batch.files) {
    message(paste("Removing the temporary batch files from directory", save.batches.dir, "with the identifier", unique.run.identifier))
    system(paste("rm ", save.batches.dir, "/", unique.run.identifier, "*", sep = ""))
  } else {
    message(paste("Keeping the temporary batch files in directory", save.batches.dir, "with the identifier", unique.run.identifier))
  }

  # Garbage collection
  gc()

  NULL

}
  



saving.hyperparameters <- function (hyper.parameters, save.batches, save.batches.dir, unique.run.identifier, verbose) {
  
  hyper.file <- NULL

  if ( save.batches ) {
    hyper.file <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-hyperparameters.RData", sep = "")
    if ( verbose ) { message(paste("Saving hyperparameters into file: ", hyper.file)) }
    save(hyper.parameters, file = hyper.file)
  }

  hyper.file

}


saving.batchlist <- function (batches, save.batches, save.batches.dir, unique.run.identifier, verbose) {

  blf <- NULL
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

  blf
}


saving.quantile.basis <- function (quantile.basis, save.batches, save.batches.dir, unique.run.identifier, verbose) {

  quantile.file <- NULL

  if (save.batches) {
    quantile.file <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-quantiles.RData", sep = "")
    if (verbose) { message(paste("Saving quantile basis into file: ", quantile.file)) }
    save(quantile.basis, file = quantile.file)
  }

  quantile.file

}

get.affinities <- function (sets, hyper.parameters, cdf, mc.cores, quantile.basis, bg.method, batches, load.batches = TRUE, save.batches, save.batches.dir, unique.run.identifier, verbose = TRUE, set.inds) {

  # load.batches = TRUE

  if ( verbose ) { message("Estimating affinities") }

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



