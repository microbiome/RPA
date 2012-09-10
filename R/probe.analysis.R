# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' probe.performance
#'
#' Provide a table of probe-level parameter estimates (affinity and stochastic noise) for RPA output.
#'
#' @param rpa.object Object of the rpa class (output from functions rpa or rpa.online)
#' @param abatch Optional: Affybatch for the rpa object, if not provided in the rpa project.
#  @param sort Sort the probes according to their stochastic noise level
#' @param sets Specify the probesets to include in the output. Default: All probesets
#' @return Data frame of probe-level parameter estimates
#' @export
#'
#' @references
#' See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # library(affydata) data(Dilution); rpa.results <- RPA.pointestimate(Dilution); tab <- probe.parameters(rpa.results); df <- df[order(abs(df$variance), decreasing = TRUE),]
#' @keywords utilities

probe.performance <- function (rpa.object, abatch = NULL, sets = NULL) {

  if (!is.null(rpa.object$abatch)) {
    # message("Using the affybatch from rpa.object")
    abatch <- rpa.object$abatch
  }

  if (is.null(sets)) {
    # Define the probesets to check
    sets <- rpa.object$sets
  }

  if (!all(sets %in% rpa.object$sets)) {
    warning("Not all sets in rpa.object, considering only the overlapping sets.")	
    sets <- intersect(sets, rpa.object$sets)
  }

  # Probe affinity effects
  af <- unlist(lapply(sets, function (set) {rpa.object[[set]]$affinity}))

  # Probe-specific noise (variance)
  s2 <- unlist(lapply(sets, function (set) {rpa.object[[set]]$tau2}))

  # PM probe indices
  pmind <- unlist(pmindex(abatch)[sets])

  # Probe effect table
  df <- data.frame(list(pmindex = pmind, affinity = af, variance = s2))

  df

}


  


