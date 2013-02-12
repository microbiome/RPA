# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' RPA.pointestimate
#'
#' Computing point estimate for the model parameters for all probe sets
#' 
#' @param abatch An AffyBatch object.
#' @param sets Specifies the probesets for which RPA estimates will be computed. Default: all probe sets.
#' @param myseed Specifies the random seed.
#' @param priors Optional list containing hyperparameters alpha and beta of the inverse Gamma prior of the probe-specific variances; alpha is a scalar, common for all probes and probeset; beta is a list where each element is a  vector corresponding to one probeset, specifying beta for each probe. Can be used to set user-specified priors for the model parameters. Not applicable for tau2.method = "var". Noninformative prior is obtained with alpha, beta approaches 0. Not used with tau2.method 'var'. Can be used to regularize the solution with small sample size, or in batch-wise online-updates with large sample size. If priors are not provided for certain probesets (NULL), default priors are used. 
#' @param epsilon Convergence tolerance. The iteration is deemed converged when the change in the d parameter is smaller than epsilon. 
#' @param cind Specifies which array in abatch is used as a reference in computing probe-level differential expression.
#' @param tau2.method Optimization method for tau2 (probe-specific variances). "robust": (default) update tau2 by posterior mean, regularized by informative priors that are identical for all probes (user-specified by setting scalar values for alpha, beta). This regularizes the solution and avoids overfitting where a single probe obtains infinite reliability. This is a potential problem in the other tau2 update methods with non-informative variance priors. The default values alpha = 2; beta = 1 are used if alpha and beta are not specified; "mode": update tau2 with posterior mean; "mean": update tau2 with posterior mean; "var": update tau2 with variance around d. Applies the fact that tau2 cost function converges to variance with large sample size
#' @param d.method Method to optimize d: "fast": (default) weighted mean over the probes, weighted by probe variances The solution converges to this with large sample size; "basic": optimization scheme to find a mode used in Lahti et al. TCBB/IEEE; relatively slow; this is the preferred method with small sample sizes.
#' @param verbose Print progress information during computation. 
#' @param bg.method Specify background correction method. Default: "rma". See bgcorrect.methods() for other options.
#' @param normalization.method Specify quantile normalization method. Default: "pmonly". See normalize.methods(Dilution) for other options.
#' @param cdf Specify an alternative CDF environment. 
#'
#' @export
#'
#' @details Calculates RPA estimates of probe reliability and differential expression between the user-specified reference array (cind) and the other arrays in the data set. The model assumes P observations for each transcript target (i.e. a probeset) with Gaussian noise which is specific for each probe (variance is specified by tau2). The mean (affinity) parameters of the Gaussian noise model cancel out in calculating probe-level differential expression.  RPA.pointestimate gives a point estimate for d and tau2. The 'prior' parameter is not applicable with tau2.method = "var". The d.method = "fast" is recommended with large sample size. tau2.method = "robust" and d.method = "fast" are recommended. With small sample size and informative priors, d.method = "basic" may be preferable, with large sample size d.method = "fast" should have considerable speedups and comparable accuracy.
#'
#' @return An instance of class 'rpa'. This is an extended list containing the following elements: d: A matrix of probesets x arrays. Specifies the estimated 'true' underlying differential gene expression signal over the arrays (vs. the reference array 'cind') for each investigated probeset. Note that the reference array is not included; tau2: A list. Each element corresponds to a probeset, and contains a vector that gives the estimated variance for each probe in that probeset. This corresponds to the parameter tau^2 in the vignette and manuscript; cind Specifies which of the arrays in the abatch (the affybatch object to be analyzed) has been used as the reference for computing probe-level differential expression; affinity: Probe affinity effects; sets: A character vector listing the investigated probesets.
#'
#' @seealso rpa.plot, rpa, set.priors, rpa2eset, RPA.preprocess, AffyBatch, rpa.fit, estimate.affinities
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords utilities methods

RPA.pointestimate <- function (abatch,
                                 sets = NULL,
                               myseed = 101, 
                               priors = NULL,
                              epsilon = 1e-2, 
                                 cind = 1,
                        tau2.method = "robust",
                             d.method = "fast",
                              verbose = TRUE,
                            bg.method = "rma",
                 normalization.method = "quantiles.robust",
                                  cdf = NULL
				 )                                      
{

# Find posterior mode for RPA model parameters d (mean) and tau2 (variances)
# and then estimate also probe affinities.

#abatch<-Dilution; sets = NULL; myseed = 101; priors = NULL; epsilon = 1e-2; cind = 1; tau2.method = "robust"; d.method = "fast"; verbose = TRUE; bg.method = "rma"; normalization.method = "quantiles.robust"; cdf = NULL

  #################################################################

  # PREPROCESSING

  #################################################################

  #Set random seed
  set.seed( myseed )

  # Set alternative CDF environment if given
  if (!is.null(cdf)) { abatch@cdfName <- cdf }
      
  # Preprocessing
  preproc <- RPA.preprocess(abatch, bg.method, normalization.method, cdf)
  
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

  ## Matrices to store the results (also including reference array)
  d.results <- array(NA, dim = c(Nsets, T))
  rownames(d.results) <- sets
  colnames(d.results) <- colnames(exprs(abatch))

  tau2.results <- vector(length = Nsets, mode = "list")  
  names(tau2.results) <- sets

  affinity.results <- vector(length = Nsets, mode = "list")  
  names(affinity.results) <- sets

  mu.real <- vector(length = Nsets, mode = "list")  
  names(mu.real) <- sets    

  # Pick the priors for this set (gives NULL if no prior has been defined)
  alpha <- priors$alpha

  for (i in 1:Nsets) {

    set <- sets[[i]]
  
    if (verbose) {message(paste("Summarizing probeset", set, ":", i, "/", Nsets, "...\n"))}

    # Find probe (pm) indices for this set
    pmindices <- preproc$set.inds[[set]]

    # Pick the priors for this set (gives NULL if no prior has been defined)
    beta  <- priors[[set]]$beta
        
    res <- rpa.fit(matrix(preproc$q[pmindices,], length(pmindices)), cind, epsilon, alpha, beta, tau2.method, d.method)
    
    #Store results
    d.results[i, ] <- res$mu # note this returns signal in original data domain
    tau2.results[[i]] <- res$tau2
    affinity.results[[i]] <- res$affinity
    mu.real[[i]] <- res$mu.real
    # mu.real can be estimated afterwards since
    # res$mu = mu.real + d and by definition d[[reference.sample]] = 0.
  }

  # Create new class instance
  # FIXME: with large data sets it exhaustive to store both abatch and preproc$q as these
  # are redundant. Even abatch is unnecessary as it is available from the input list already.
  rpa.res <- new("rpa", list(d = d.results, mu.real = mu.real, tau2 = tau2.results, affinity = affinity.results, cind = cind, sets = sets, data = preproc$q, cdf = cdf, abatch = abatch))

  # return result object
  rpa.res
}

