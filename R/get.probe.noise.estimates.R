# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. 
# All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' get.probe.noise.estimates
#' Fetch probe-level noise estimates from an rpa object
#' 
#' @param rpa.res An rpa object.
#' @param sets Probesets to check. 
#' @param normalization Optional normalization for probe noise estimates. The higher the value, the higher the probe-level noise. By default, probe-level variances of the RPA model are returned. Other options include:
#'    "withinset.weights": The relative weight of a probe within probeset
#'       is determined by the relative noise of the probe with respect to
#'       the other probes in the same probeset. This option returns the
#'       inverse of probe-specific weights within each probeset. This can
#'       be used to normalize probe-level weights to improve comparability
#'       across probesets.
#'       "withinset.relative": The detected probe-level noise can be
#'         coupled with overall signal levels of the probeset. This option
#'         provides an estimate of probe-wise standard deviation
#'         normalized by the standard deviation of the probeset-level
#'         signal d.
#'	 "withinset.categorical": In some applications it can be
#'    	   sufficient to investigate the relative order of the probes,
#'    	   ignoring the parameter estimates. This option indexes the
#'    	   probes according to their reliability within each
#'    	   probeset. Probes with higher indices are more noisy.  
#' @param verbose Print progress information during computation.
#'
#' @details The normalization options are included to improve comparability across probesets. The higher the variance, the more noisy the probe. Inverse of the variance, can be used to quantitate probe reliability. Note that the relative weight of a probe within probeset is determined by the relative noise of the probe with respect to the other probes in the same probeset. Comparison of probe-specific variances across probesets may benefit from normalization of this effect. Therefore optional normalizations for probe noise estimation are provided.
#'
#' @return A list. Each element corresponds to one probeset (of the input object). The element lists noise estimates for each probe within the probeset.
#'
#' @seealso RPA.pointestimate
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples # 
#' @keywords utilities
#'

get.probe.noise.estimates <- function (rpa.res, sets = NULL, normalization = NULL, verbose = FALSE) {

  if (is.null(sets)) {sets <- names(rpa.res)}
  
  if ( is.null( normalization ) ) {
    # No normalization: return probe variances
    rpa.res$tau2 <- rpa.res$tau2

  } else if (normalization == "withinset.categorical") {

    for (set in sets) {
      # Order probes within the probeset by their noise level
      # smallest index means least noise
      if (verbose) {cat(set)}
      o <- order(rpa.res$tau2[[set]])
      rpa.res$tau2[[set]][o] <- 1:length(o)      
    }
    
  } else if (normalization == "withinset.relative") {

    for (set in sets) {
      # compare probe standard deviation to probeset signal standard
      # deviation
      # within the probeset
      # smallest means least noise
      if (verbose) {cat(set)}
      vars <- rpa.res$tau2[[set]]
      rpa.res$tau2[[set]] <- sqrt(vars)/sqrt(var(rpa.res$d[set,]))
    }
  } else if (normalization == "withinset.weights") {

    for (set in sets) {
      # inverse of the weights used to compute
      # probeset-level summary signal d
      # within the probeset
      # smallest means least noise (largest weight = smallest inverse weight)
      if (verbose) {cat(set)}
      vars <- rpa.res$tau2[[set]]
      rpa.res$tau2[[set]] <- 1/((1 / vars) / sum(1 / vars))
    }
  }

  rpa.res$tau2[sets]

}
