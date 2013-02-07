# This file is a part of the RPA (Robust Probabilistic Averaging)
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' rpa.summarize
#'
#' @param dat Original data: probes x samples.
#' @param affinities Probe affinities
#' @param variances Probe variances
#'
#' @details Summarizes the probes in a probe set according to the RPA model based on the given affinity and variance parameters.
#'
#' @returns A vector. Probeset-level summary signal.
#'
#' @seealso rpa, RPA.pointestimate
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # res <- rpa.summarize(dat, affinities, variances)
#' @keywords utilities

rpa.summarize <- function (dat, affinities, variances) {

  # FIXME: add hyperparameter estimation for the case where 
  # affinities & variances are NULL

  # Impute if there are missing values
  dat <- rpa.impute(dat)
  if (is.null(colnames(dat))) {colnames(dat) <- 1:ncol(dat)}

  # Accommodate single-probe probesets
  if (nrow(dat) == 1) {  

    mu <- as.vector(dat)
    mu <- mu - affinities
    names(mu) <- colnames(dat)

  } else {

    # Remove affinities from raw signal before summarization

    # Since probe-specific variance is now known (from the estimation above), 
    # the probeset-level signal
    # estimate is obtained as a weighted sum of the 
    # probes, weighted by the probe-specific variances

    #mu <- d.update.fast(dat - affinities, variances)
    mu <- d.update.fast(dat, variances) # ignore affinities so far
  }
    
  mu

}
