# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' get.probe.matrix
#'
#'
#' @param cels List of CEL files to preprocess
#' @param cdf Specify an alternative CDF environment
#' @param quantile.basis Pre-calculated basis for quantile normalization
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
  
  if (!is.null(batch)) {

      # Assuming that the bg correction + quantile normalization have
      # been already calculated for quantile.basis, which is here
      # simply allocated for each array
    
      if (verbose) { message("Set quantile data on each array") }

      q <- apply(batch, 2, function (o) { quantile.basis[o] }) 

      if (verbose) { message("...Done.") }
      
  } else {

      # background-corrected, quantile normalized log2 probes x samples matrix
      q <- RPA.preprocess(abatch = NULL, bg.method = bg.method, normalization.method = normalization.method, cdf = cdf, cel.files = cels, quantile.basis = quantile.basis)$q	

  }

  q

}

