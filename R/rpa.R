# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' rpa
#' RPA for preprocessing.
#'
#' @param abatch An AffyBatch object.
#' @param sets Probesets for which RPA will be computed. 
#' @param priors An 'rpa.priors' object. Can be used to set user-specified priors for the model parameters. Not used tau2.method = "var". The prior parameters alpha and beta are prior parameters for inverse Gamma distribution of probe-specific variances. Noninformative prior is obtained with alpha, beta -> 0.  Not used with tau2.method 'var'. Scalar alpha and beta specify an identical inverse Gamma prior for all probes, which regularizes the solution. Can be also specified as lists, each element corresponding to one probeset.
#' @param epsilon Convergence tolerance. The iteration is deemed converged when the change in all parameters is < epsilon.
#' @param cind Specify reference array for computing probe-level differential expression. Default: cind = 1. Note that if exclude.reference.array = TRUE the expression value for the reference array (cind) will be excluded in the output. Note that all values of the reference array are 0 since they indicate the differential expression of the reference array against itself.
#' @param tau2.method Optimization method for tau2 (probe-specific variances). This parameter is denoted by tau^2 in the vignette and manuscript
#'
#'	"robust": (default) update tau2 by posterior mean,
#'		regularized by informative priors that are identical
#'		for all probes (user-specified by
#'		setting scalar values for alpha, beta). This
#'		regularizes the solution, and avoids overfitting where
#'		a single probe obtains infinite reliability. This is a
#'	        potential problem in the other tau2 update
#'	        methods with non-informative variance priors. The
#'		default values alpha = 2; beta = 1 are
#'	        used if alpha and beta are not specified.
#'   
#'      "mode": update tau2 with posterior mean
#'
#'	"mean": update tau2 with posterior mean
#'	
#'	"var": update tau2 with variance around d. Applies the fact
#'             that tau2 cost function converges to variance with
#'               large sample sizes. 
#'
#' @param d.method Method to optimize d.
#'
#'        "fast": (default) weighted mean over the probes, weighted by
#'		probe variances The solution converges to this with
#'		large sample size.
#'
#'      "basic": optimization scheme to find a mode used in Lahti et
#'        	 al. TCBB/IEEE; relatively slow; this is the preferred 
#'		 method with small sample sizes.
#'     
#' @param verbose Print progress information during computation.
#' @param bg.method Specify background correction method. Default: "rma". See bgcorrect.methods() for other options.
#' @param normalization.method Specify quantile normalization method. Default: "pmonly". See normalize.methods(Dilution) for other options.
#' @param cdf Specify an alternative CDF environment. Default: none.
#' @param cel.files List of CEL files to preprocess.
#' @param cel.path Path to CEL file directory.
#' @param probe.parameters A list, each element corresponding to a probe set. Each probeset element has the following elements: affinity, variance. Each of these two elements contains a vector over the probeset probes, specifying the affinities and variances for the probes according to the RPA model. If given, overrides the priors.
#' @param mc.cores Number of cores for parallelized processing. 
#'
#' @details Returns an expressionSet object preprocessed with RPA. If 'cind' is not specified, uses the first array of affybatch as the reference. RPA preprocessing function. Gives an estimate of the probeset-level mean parameter d of the RPA model, and returns these in an expressionSet object. The choices tau2.method = "robust" and d.method = "fast" are recommended. With small sample size and informative prior, d.method = "basic" may be preferable. For very large expression data collections, see rpa.online function.
#'
#' @return An instance of the 'expressionSet' class
#'
#' @seealso rpa.online, RPA.pointestimate, set.priors, AffyBatch, ExpressionSet, estimate.affinities, rpa.fit
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # eset <- rpa(abatch)
#' @keywords methods

rpa <- function (abatch = NULL,
                    sets = NULL,
                  priors = NULL,
                 epsilon = 1e-2, 
                    cind = 1,
           tau2.method = "robust",
                d.method = "fast",
                 verbose = FALSE,
               bg.method = "rma",
    normalization.method = "quantiles.robust",
                     cdf = NULL,
               cel.files = NULL, 
                cel.path = NULL, 
	probe.parameters = NULL, 
	        mc.cores = 1) 
{

  # abatch = NULL; sets = NULL; priors = NULL; epsilon = 1e-2; cind = 1; tau2.method = "robust"; d.method = "fast"; verbose = FALSE; bg.method = "rma"; normalization.method = "quantiles.robust"; cdf = NULL; cel.files = NULL; cel.path = "CEL"; probe.parameters = NULL; mc.cores = 1

  # Background correction, normalization and logging
  preproc <- RPA.preprocess(abatch, bg.method, normalization.method, cdf = cdf, cel.files = cel.files, cel.path = cel.path)	

  if ( verbose ) { message("Calculating Expression") }

  # Summarization including hyperparameter estimation, possibly with priors.
  mat <- summarize.batch(preproc$q, preproc$set.inds, 
    	   	         priors, probe.parameters, 
			 cind, epsilon,			
			 verbose, mc.cores) 
  
  # Coerce expression values in the rpa object into an ExpressionSet object
  # and return expression set
  new("ExpressionSet", assayData = list(exprs = mat)) 
  
}



