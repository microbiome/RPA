#' @title get.probe.matrix
#' @description Get probe matrix.
#'
#' @param cels List of CEL files to preprocess
#' @param cdf Specify an alternative CDF environment
#' @param quantile.basis Pre-calculated basis for quantile normalization in log2 domain
#' @param bg.method Specify background correction method. See bgcorrect.methods() for options.
#' @param normalization.method normalization method
#' @param batch batch
#' @param verbose Print progress information during computation
#'
#' @details Returns background-corrected, quantile normalized log2 probes x samples matrix
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords methods

get.probe.matrix <- function (cels, cdf = NULL, quantile.basis, bg.method = "rma", normalization.method = "quantiles", batch = NULL, verbose = TRUE) {
  
  # cels = batch.cels; quantile.basis = probe.parameters$quantile.basis; 

  if (!is.null(batch)) {

      # Assuming that the bg correction + quantile normalization have
      # been already calculated for quantile.basis, which is here
      # simply allocated for each array
    
      if (verbose) { message("Set quantile data on each array") }

      q <- apply(batch, 2, function (o) { quantile.basis[o] }) 

      if (verbose) { message("...Done.") }
      
  } else {

      # background-corrected, quantile normalized log2 probes x samples matrix
      # NOTE: this requires quantile.basis in original domain, not log2 !!!
      q <- RPA.preprocess(abatch = NULL, bg.method = bg.method, normalization.method = normalization.method, cdf = cdf, cel.files = cels, quantile.basis = 2^quantile.basis)$q	

  }

  q

}

