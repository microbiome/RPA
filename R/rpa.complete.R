#' @title Complete RPA preprocessing
#' @description RPA preprocessing, also returns probe parameters.
#' @param abatch An AffyBatch object.
#' @param sets Probesets for which RPA will be computed. 
#' @param epsilon Convergence tolerance. The iteration is deemed converged when the change in all parameters is < epsilon.
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
#' @param probe.parameters A list, each element corresponding to a probe set. Each probeset element has the following optional elements: affinity (affinity), tau2 (variance), alpha (shape prior), betas (scale prior). Each of these elements contains a vector over the probeset probes, specifying the probe parameters according to the RPA model. If variance is given, it overrides the priors. Can be also used to set user-specified priors for the model parameters. Not used tau2.method = "var". The prior parameters alpha and beta are prior parameters for inverse Gamma distribution of probe-specific variances. Noninformative prior is obtained with alpha, beta -> 0.  Not used with tau2.method 'var'. Scalar alpha and beta specify an identical inverse Gamma prior for all probes, which regularizes the solution. Can be also specified as lists, each element corresponding to one probeset. Can also include quantile.basis 
#' @param mc.cores Number of cores for parallelized processing. 
#' @param summarize.with.affinities Use affinity estimates in probe summarization step. Default: FALSE.
#'
#' @details RPA preprocessing function. Gives an estimate of the probeset-level mean parameter d of the RPA model, and returns these in an expressionSet object. The choices tau2.method = "robust" and d.method = "fast" are recommended. With small sample size and informative prior, d.method = "basic" may be preferable. For very large expression data collections, see rpa.online function.
#'
#' @return List with preprocessed expression matrix, corresponding probe parameters, AffyBatch and CDF
#'
#' @export
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # eset <- rpa(abatch)
#' @keywords methods

rpa.complete <- function (abatch = NULL,
                    sets = NULL,
                 epsilon = 1e-2, 
           tau2.method = "robust",
                d.method = "fast",
                 verbose = FALSE,
               bg.method = "rma",
    normalization.method = "quantiles.robust",
                     cdf = NULL,
               cel.files = NULL, 
                cel.path = NULL, 
	probe.parameters = list(), 
	        mc.cores = 1, 
	  summarize.with.affinities = FALSE
		) 
{


   # Convert probe parameters to list if needed:
  probe.parameters <- probe.parameters.tolist(probe.parameters)

  # Background correction, normalization and logging
  preproc <- RPA.preprocess(abatch = abatch, 
  	     		    bg.method = bg.method, 
			    normalization.method = normalization.method, 
			    cdf = cdf, 
			    cel.files = cel.files, 
			    cel.path = cel.path, 
		      quantile.basis = probe.parameters$quantile.basis)	

  # Store the applied quantile.basis		      
  message("Setting quantile basis")
  probe.parameters$quantile.basis <- preproc$quantile.basis

  # Summarization including hyperparameter estimation, possibly with priors.
  res <- summarize.batch(      q = preproc$q, 
      	 	        set.inds = preproc$set.inds, 
	        probe.parameters = probe.parameters, 
			 epsilon = epsilon,			
			 verbose = verbose, 
			mc.cores = mc.cores, 
       summarize.with.affinities = summarize.with.affinities) 

  # Convert expression matrix into expressionSet format before returning
  list(eset = new("ExpressionSet", assayData = list(exprs = res$exprs)), probe.parameters = res$probe.parameters, abatch = abatch, cdf = cdf, probedata = preproc$q)
		
}



