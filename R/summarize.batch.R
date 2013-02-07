# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' summarize.batch
#'
#' @param q Background corrected, quantile-normalized, log2 probes x samples matrix
#' @param set.inds Indices for each probeset, corresponding to q matrix
#' @param priors Optional probe-specific priors
#' @param probe.parameters A list, each element corresponding to a probe set. Each probeset element has the following elements: affinity, variance. Each of these two elements contains a vector over the probeset probes, specifying the affinities and variances for the probes according to the RPA model. If given, overrides the priors.
#' @param cind Specify reference array for computing probe-level differential expression. Default: cind = 1. Note that if exclude.reference.array = TRUE the expression value for the reference array (cind) will be excluded in the output. Note that all values of the reference array are 0 since they indicate the differential expression of the reference array against itself.
#' @param epsilon Convergence tolerance. The iteration is deemed converged when the change in all parameters is < epsilon.
#' @param verbose Print progress information during computation.
#' @param mc.cores Number of cores for parallel processing
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords methods

summarize.batch <- function (q, set.inds, priors = NULL, probe.parameters = NULL, cind = 1, epsilon, verbose = FALSE, mc.cores = 1) {

  # q <- preproc$q; set.inds <- preproc$set.inds

  sets <- names(set.inds)
  Nsets <- length( sets ) 

  # Pick each probe set separately    
  q <- mclapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) }, mc.cores = mc.cores) 
  names(q) <- sets

  if (!is.null(probe.parameters)) {
 
    # Get summary estimates for each probeset with posterior variances
    # (>25-fold speedup with sapply!)
    if ( verbose ) { message("Summarizing...") }       
    emat <- t(sapply(mclapply(sets, function(set) {    	   
    	       rpa.summarize(q[[set]], 
	       		     probe.parameters$affinity[[set]],
			     probe.parameters$variance[[set]])
	       }, mc.cores = mc.cores), identity))

  } else {

     # Estimate model parameters from the data
     emat <- t(sapply(mclapply(sets, function(set) {
	       rpa.fit(dat = q[[set]], 
	       	       alpha = priors[[set]]$alpha, 
		       beta = priors[[set]]$beta)$mu
	     }, mc.cores = mc.cores), identity))
  }

  emat
  
}


