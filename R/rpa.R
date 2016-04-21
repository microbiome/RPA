#' @title rpa
#' @description Wrapper for RPA preprocessing.
#' @param abatch An AffyBatch object.
#' @param verbose Print progress information during computation.
#' @param bg.method Specify background correction method. Default: "rma". See bgcorrect.methods() for other options.
#' @param normalization.method Specify quantile normalization method. Default: "pmonly". See normalize.methods(Dilution) for other options.
#' @param cdf Specify an alternative CDF environment. Default: none.
#' @param cel.files List of CEL files to preprocess.
#' @param cel.path Path to CEL file directory.
#' @param probe.parameters A list, each element corresponding to a probe set. Each probeset element has the following optional elements: mu (affinity), tau2 (variance), alpha (shape prior), beta (scale prior). Each of these elements contains a vector over the probeset probes, specifying the probe parameters according to the RPA model. If variance is given, it overrides the priors. Can be also used to set user-specified priors for the model parameters. Not used tau2.method = "var". The prior parameters alpha and beta are prior parameters for inverse Gamma distribution of probe-specific variances. Noninformative prior is obtained with alpha, beta -> 0.  Not used with tau2.method 'var'. Scalar alpha and beta specify an identical inverse Gamma prior for all probes, which regularizes the solution. Can be also specified as lists, each element corresponding to one probeset. May also include quantile.basis
#' @param mc.cores Number of cores for parallelized processing. 
#' @param summarize.with.affinities Use affinity estimates in probe summarization step. Default: FALSE.
#'
#' @details RPA preprocessing function. Gives an estimate of the probeset-level mean parameter d of the RPA model, and returns these in an expressionSet object. The choices tau2.method = "robust" and d.method = "fast" are recommended. With small sample size and informative prior, d.method = "basic" may be preferable. For very large expression data collections, see rpa.online function.
#'
#' @return Preprocessed expression matrix in expressionSet format
#'
#' @seealso rpa.online, AffyBatch, ExpressionSet, estimate.affinities, rpa.fit
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # eset <- rpa(abatch)
#' @keywords methods

rpa <- function (abatch = NULL,
                 verbose = FALSE,
               bg.method = "rma",
    normalization.method = "quantiles.robust",
                     cdf = NULL,
               cel.files = NULL, 
                cel.path = NULL, 
	probe.parameters = NULL, 
	        mc.cores = 1, 
	  summarize.with.affinities = FALSE
		) 
{

  # bg.method = "none"; cdf = CDF_NAME; verbose = FALSE; normalization.method = "quantiles.robust"; cel.files = NULL; cel.path = NULL; probe.parameters = NULL; mc.cores = 1; summarize.with.affinities = FALSE

  if (ncol(abatch) <= 1) {
    warning("RPA is a multi-array preprocessing method. The input affybatch has at most a single array, however. Returning expressionSet with no background correction or normalization.")
    eset <- affy::rma(abatch, background = FALSE, normalize = FALSE)
    return(eset)
  }

  # RPA analysis	
  res <- rpa.complete(abatch = abatch,
                 verbose = verbose,
               bg.method = bg.method,
    normalization.method = normalization.method,
                     cdf = cdf,
               cel.files = cel.files, 
                cel.path = cel.path, 
	probe.parameters = probe.parameters, 
	        mc.cores = mc.cores, 
	  summarize.with.affinities = summarize.with.affinities
		) 

  # Return ExpressionSet object		
  res$eset
  
}


