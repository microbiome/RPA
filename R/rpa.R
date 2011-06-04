# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  
rpa <- function (abatch = NULL,
                    sets = NULL,
                  priors = NULL,
                 epsilon = 1e-2, 
                    cind = 1,
           sigma2.method = "robust",
                d.method = "fast",
                 verbose = FALSE,
               bg.method = "rma",
    normalization.method = "quantiles.robust",
                     cdf = NULL,
               cel.files = NULL, 
                cel.path = NULL) 
{

  # RPA preprocessing wrapper

  # Background correction, normalization and logging
  preproc <- RPA.preprocess(abatch, bg.method, normalization.method, cdf = cdf, cel.files = cel.files, cel.path = cel.path)	

  message("Calculating Expression")
  
  # Check names and number for the investigated probesets
  # if my.sets not specified, take all sets in abatch
  if (is.null(sets)) { sets <- names(preproc$set.inds) } 
  Nsets <- length( sets ) 

  ## Matrices to store results
  mu.results <- array(NA, dim = c(Nsets, ncol( preproc$q )))
  rownames(mu.results) <- sets
  colnames(mu.results) <- colnames(preproc$q)

  for (i in 1:Nsets) {

    set <- sets[[i]]
  
    if (verbose) { message(paste("Summarizing probeset", set, ":", i, "/", Nsets, "...\n")) }

    # Find probe (pm) indices for this set
    pmindices <- preproc$set.inds[[set]] #pmindex(abatch, set)[[1]]
      
    # Pick the priors for this set
    alpha <- priors[[set]]$alpha 
    beta  <- priors[[set]]$beta
    
    # Calculate RPA 
    # (note: in the new version this is done in the original data domain, 
    # including the reference array)
    dat <- matrix(preproc$q[pmindices,], length(pmindices))
    mu.results[i, ] <- rpa.fit(dat, cind, epsilon, alpha, beta, sigma2.method, d.method)$mu
    # Store the results (only mean parameter and in the original data domain)    
     
  }

  # Coerce expression values in the rpa object into an ExpressionSet object
  # and return expression set
  new("ExpressionSet", assayData = list(exprs = mu.results)) 
  
}

