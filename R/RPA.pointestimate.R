# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


RPA.pointestimate <- function (abatch,
                                 sets = NULL,
                               myseed = 101, 
                               priors = NULL,
                              epsilon = 1e-2, 
                                 cind = 1,
                        sigma2.method = "robust",
                             d.method = "fast",
                              verbose = TRUE,
                            bg.method = "rma",
                 normalization.method = "quantiles.robust",
                                  cdf = NULL,
                                alpha = NULL,
                                 beta = NULL,
		      quantile.n = 50
				 )                                      
{


# Find posterior mode for RPA model parameters d (mean) and sigma2 (variances)
# and then estimate also probe affinities.

  #################################################################

  # PREPROCESSING

  #################################################################

  #Set random seed
  set.seed( myseed )

  # Set alternative CDF environment if given
  if (!is.null(cdf)) { abatch@cdfName <- cdf }
      
  # Preprocessing
  preproc <- RPA.preprocess(abatch, cind, bg.method, normalization.method, cdf, quantile.n = quantile.n)
  
  #################################################################

  # ESTIMATE PROBE RELIABILITY AND DIFFERENTIAL GENE EXPRESSION

  #################################################################

  # Number of arrays 
  T <- ncol(exprs(abatch))

  # Check names and number for the investigated probesets
  # if my.sets not specified, take all sets in abatch
  if ( is.null(sets) ) { sets <- geneNames(abatch) } 

  if (!all(sets %in% probeNames(abatch))) {
    warning("some probesets (sets) not available in abatch!")
    sets <- sets[sets %in% probeNames(abatch)]
  } else {}
  
  Nsets <- length(sets) 

  ## Matrices to store the results
  d.results <- array(NA, dim = c(Nsets, T))
  rownames(d.results) <- sets
  colnames(d.results) <- colnames(exprs(abatch))

  sigma2.results <- vector(length = Nsets, mode = "list")  
  names(sigma2.results) <- sets

  affinity.results <- vector(length = Nsets, mode = "list")  
  names(affinity.results) <- sets

  mu.real <- vector(length = Nsets, mode = "list")  
  names(mu.real) <- sets    

  if (!is.null(priors) && (!is.null(alpha) || !is.null(beta))) {
    alpha <- beta <- NULL
    warning("priors parameter is overriding alpha, beta when both are provided in the function input")
  }

  for (i in 1:Nsets) {

    set <- sets[[i]]
  
    if (verbose) {message(paste("Summarizing probeset", set, ":", i, "/", Nsets, "...\n"))}

    # Find probe (pm) indices for this set
    pmindices <- preproc$set.inds[[set]]
    
    # Number of probes in this probeset
    P <- length(pmindices)

    # Pick the priors for this set (gives NULL if no prior has been defined)
    if (!is.null(priors)) {
      alpha <- priors[[set]]$alpha 
      beta  <- priors[[set]]$beta
    }

    res <- rpa.fit(preproc$q[pmindices,], cind, epsilon, alpha, beta, sigma2.method, d.method)
    
    #Store results
    d.results[i, ] <- res$mu # note this returns signal in original data domain
    sigma2.results[[i]] <- res$sigma2
    affinity.results[[i]] <- res$affinity
    mu.real[[i]] <- res$mu.real
    # mu.real can be estimated afterwards since
    # res$mu = mu.real + d and by definition d[[reference.sample]] = 0.
  }

  # Create new class instance
  # FIXME: with large data sets it exhaustive to store both abatch and preproc$q as these
  # are redundant. Even abatch is unnecessary as it is available from the input list already.
  rpa.res <- new("rpa", list(d = d.results, mu.real = mu.real, sigma2 = sigma2.results, affinity = affinity.results, cind = cind, sets = sets, data = preproc$q, cdf = cdf, abatch = abatch))

  # return result object
  rpa.res
}

