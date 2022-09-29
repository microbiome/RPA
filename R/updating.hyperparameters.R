#' @title updating hyperparameters
#' @description Hyperparameter update.
#'
#' @param q probes x samples matrix
#' @param set.inds Probe set indices
#' @param verbose Print progress information
#' @param mc.cores Number of cores for parallel computation
#' @param alpha alpha hyperparameter 
#' @param betas beta hyperparameters
#' @param epsilon Convergence parameter
#'
#' @return List with the following elements: alpha, betas, s2s (variances) 
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples # 
#' @keywords utilities

updating.hyperparameters <- function (q, set.inds, verbose, mc.cores = 1, alpha, betas, epsilon) {

    # Get probes x samples matrix of probe-wise fold-changes
    # Select one of the arrays as a reference at
    # random for each batch. Choice of the reference array does not notably affect
    # the results in experiments as the control effect is marginalized
    # out in the treatment. Note that in rpa.online implementation,
    # cind is specific to each batch but it is only used to in
    # hyperparameter estimation step to cancel probe affinity effects;
    # in probeset summarization no reference sample is needed. Whether
    # cind is the same for the overall data collection or
    # batch-specific should not notably affect the results, either.

    cind <- sample(ncol(q), 1)
    q <- matrix(q[, -cind] - q[, cind], nrow(q))

    if ( verbose ) { message("Get probes x samples matrices for each probeset") }
    sets <- names(set.inds)

    q <- mclapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) }, mc.cores = mc.cores)
    names(q) <- sets

    # Get variance point estimate. Then beta can be solved, as alpha is given based on batch size
    if ( verbose ) { message("Update probe parameters") }

    pars <- mclapply(sets, function (set) { 
    	 
	 # For a single probeset
         estimated <- RPA.iteration(t(q[[set]]), epsilon, alpha, betas[[set]]);
    	 
         return(list(alpha = estimated$alpha,
	     	      beta = estimated$beta,
    	 	       s2s = estimated$tau2))

    }, mc.cores = mc.cores)

    # Alpha is same for all probesets
    alpha <- pars[[1]]$alpha

    # Pick betas and variances into their own vectors
    betas <- mclapply(pars, function (x) {x$beta}, mc.cores = mc.cores)
    names(betas) <- sets

    s2s   <- mclapply(pars, function (x) {x$s2s}, mc.cores = mc.cores)    
    names(s2s) <- sets

    list(alpha = alpha, betas = betas, s2s = s2s)

}
