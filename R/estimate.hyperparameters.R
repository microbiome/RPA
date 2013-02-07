# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2010-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' estimate.hyperparameters
#' Hyperparameter estimation
#'
#' @param sets Probesets to handle. All probesets by default.
#' @param priors User-defined priors
#' @param batches Data batches for online learning
#' @param cdf CDF probeset definition file
#' @param quantile.basis Basis for quantile normalization
#' @param bg.method Background correction method
#' @param epsilon Convergence parameter
#' @param load.batches Logical. Load preprocessed data whose identifiers are picked from names(batches). Assuming that the same batch list (batches) was used to create the files in online.quantiles function.
#' @param save.hyperparameter.batches Save hyperparameters for each batch into files using the identifiers with batch name with -hyper.RData suffix.
#' @param mc.cores Number of cores for parallel computation
#' @param verbose Print progress information
#' @param normalization.method Normalization method
#' @param save.batches.dir Specify the output directory for temporary batch saves.
#' @param unique.run.identifier Define identifier for this run for naming the temporary batch files. By default, a random id is generated.
#' @param set.inds Probeset indices
#'
#'@return alpha: Hyperparameter alpha (same for all probesets); betas: Hyperparameter beta (probe-specific); variances: Probe-specific variances (beta/alpha)
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples # 
#' @keywords utilities

estimate.hyperparameters <- function (sets = NULL, 
			            priors = list(alpha = 2, beta = 1), 
			    	      batches, cdf = NULL, quantile.basis, 
				      bg.method = "rma", 
				      epsilon = 1e-2,
                                      load.batches = FALSE,
                                      save.hyperparameter.batches = FALSE, 
				      mc.cores = 1, 
				      verbose = TRUE, 
				      normalization.method = "quantiles",
				      save.batches.dir = ".", 
				      unique.run.identifier = NULL, set.inds)
{

  # Hyperparameter estimation through batches  

  if ( is.null(sets) ) { sets <- names(set.inds) }
  
  # Initialize hyperparameters
  # Note: alpha is scalar and same for all probesets 
  # alpha <- alpha + N/2 at each batch
  if (verbose) { message("Initialize priors") }
  alpha <- priors$alpha # initialize 
  betas <- mclapply(sets, function (set) { 
    rep(priors$beta, length(set.inds[[set]])) 
  }, mc.cores = mc.cores)
  names(betas) <- sets

  for (i in 1:length(batches)) {

    if (verbose) {
      message(paste("Updating hyperparameters; batch", i, "/", length(batches)))
    }
    
    # Load batch with bgc, ordered data
    batch <- NULL
    if (load.batches) {
      batch.file <- paste(save.batches.dir, "/", unique.run.identifier, "-", names(batches)[[i]], ".RData", sep = "")      
      if (verbose) { message(paste("Load batch from file:", batch.file)) }
      load(batch.file) # batch
    } 

    # Get background corrected, quantile normalized, log2 probe-level matrix
    if ( verbose ) { message("Pick probe-level values") }    
    q <- get.probe.matrix(cels = batches[[i]], cdf, quantile.basis,
                          bg.method, normalization.method, batch,
                          verbose = verbose)

    # --------------------------------------------------------------------------------------------------
    
    # Update hyperparameters
    hp <- updating.hyperparameters(q, set.inds, verbose, mc.cores = 1, alpha, betas, epsilon)
    alpha <- hp$alpha
    betas <- hp$betas
    s2s   <- hp$s2s

    bf <- saving.hyperparameter.batches(alpha, betas, save.hyperparameter.batches, save.batches.dir, unique.run.identifier, names(batches)[[i]], verbose)
        
  }

  # Get final estimated variances for each probeset based on hyperparameter posteriors
  # variances <- mclapply(betas, function (beta) {beta/alpha}, mc.cores = mc.cores)

  list(alpha = alpha, betas = betas, tau2 = s2s)  

}


saving.hyperparameter.batches <- function (alpha, betas, save.hyperparameter.batches, save.batches.dir, unique.run.identifier, nam, verbose) {

    if (save.hyperparameter.batches) {
      batch.file <- paste(save.batches.dir, "/", unique.run.identifier, nam, "-hyper.RData", sep = "")
      if ( verbose ) { message(paste("Save hyperparameters into file:", batch.file)) }
      save(alpha, betas, file = batch.file)
    }
}


