#' @title frpa
#' @description Frozen-RPA preprocessing using precalculated probe parameters.
#' @param abatch An AffyBatch object.
#' @param probe.parameters A list with tau2 (probe variance), quantile.basis (basis for quantile normalization in log2 domain), and optionally affinity (probe affinities). The probe.parameters$tau2 and probe.parameters$affinity are lists, each element corresponding to a probeset and containing a parameter vector over the probes. The quantile.basis is a vector over the probes, the probes need to be listed in the same order as in tau2 and affinity. probe.parameters can be optionally provided as a data frame.
#' @param verbose Print progress information during computation.
#' @param cdf Specify an alternative CDF environment. Default: none.
#' @param cel.files List of CEL files to preprocess.
#' @param cel.path Path to CEL file directory.
#' @param mc.cores Number of cores for parallelized processing. 
#' @param summarize.with.affinities Use affinity estimates in probe summarization step. Default: FALSE.
#'
#' @details fRPA function to preprocess Affymetrix CEL files with RPA using precalculated (frozen) probe parameters.
#'
#' @return Preprocessed expression matrix in expressionSet format
#'
#' @seealso rpa, AffyBatch, ExpressionSet
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # eset <- frpa(abatch, probe.parameters)
#' @keywords methods

frpa <- function (abatch = NULL,
	probe.parameters = NULL,
                 verbose = FALSE,
                     cdf = NULL,
               cel.files = NULL, 
                cel.path = NULL, 
	        mc.cores = 1, 
	  summarize.with.affinities = FALSE
		) 
{

  # RPA analysis	
  eset <- rpa(abatch = abatch,
                 verbose = verbose,
                     cdf = cdf,
               cel.files = cel.files, 
                cel.path = cel.path, 
	probe.parameters = probe.parameters, 
	        mc.cores = mc.cores, 
	  summarize.with.affinities = summarize.with.affinities
		) 

  # Return ExpressionSet object		
  eset
  
}



