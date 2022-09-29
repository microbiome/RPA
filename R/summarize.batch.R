#' @title summarize.batch
#' @description Summarize batch.
#' @param q Background corrected, quantile-normalized, log2 probes x samples matrix
#' @param set.inds Indices for each probeset, corresponding to q matrix
#' @param probe.parameters A list, each element corresponding to a probe set. Each probeset element has the following elements: affinity, variance and optionally alpha and beta priors. Each of these elements contains a vector over the probeset probes, specifying the probe parameters according to the RPA model. If variances are given, that overrides the priors.
#' @param epsilon Convergence tolerance. The iteration is deemed converged when the change in all parameters is < epsilon.
#' @param verbose Print progress information during computation.
#' @param mc.cores Number of cores for parallel processing
#' @param summarize.with.affinities Use affinity estimates in probe summarization step. Default: FALSE.
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords methods

summarize.batch <- function (q, set.inds, probe.parameters = list(), epsilon, verbose = FALSE, mc.cores = 1, summarize.with.affinities = FALSE) {

  # q = preproc$q; set.inds = preproc$set.inds; mc.cores = 1

  if (verbose) {message("Summarizing batch")}

  sets <- names(set.inds)
  Nsets <- length( sets ) 
  sample.names <- colnames(q)

  # Pick each probe set separately    
  q <- mclapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) }, mc.cores = mc.cores) 
  names(q) <- sets

  if (!is.null(probe.parameters$tau2)) {
 
    # Get summary estimates for each probeset using pre-calculated
    # posterior (affinities) and variances (>25-fold speedup with sapply!)
    if ( verbose ) { message("Summarizing...") }       
    emat <- t(sapply(mclapply(sets, function(set) {    	   
    	       rpa.summarize(q[[set]], 
	       		     probe.parameters$affinity[[set]],
			     probe.parameters$tau2[[set]])
	       }, mc.cores = mc.cores), identity))

  } else {

     # Estimate model parameters from the data and retrieve summaries
     if ( verbose ) { message("Estimating parameters and summarizing...") }       
     res <- mclapply(sets, function(set) {
	       rpa.fit(dat = q[[set]], 
	       	     alpha = probe.parameters$alpha, 
		      beta = probe.parameters$beta[[set]],  
    summarize.with.affinities = summarize.with.affinities)
	     }, mc.cores = mc.cores)

      emat <- t(sapply(res, function (x) {x$mu}))

      # Add estimated parameters to probe.parameters
      tau2 <- lapply(res, function (x) {x$tau2})
      names(tau2) <- sets
      probe.parameters$tau2 <- tau2

      affinity <- lapply(res, function (x) {x$affinity})
      names(affinity) <- sets
      probe.parameters$affinity <- affinity

  }

  rownames(emat) <- sets
  colnames(emat) <- sample.names
  
  list(exprs = emat, probe.parameters = probe.parameters)
  
}


