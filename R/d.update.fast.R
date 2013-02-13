# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' d.update.fast
#'
#' Description: computes weighted average over the probes, 
#' weighted by their inverse probe-specific variances.
#'
#' @param St probes x samples data matrix
#' @param s2 variances for the probes
#'
#' @details Returns summarized probeset-level weighted average
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords methods

d.update.fast <- function (St, s2) {

  # With large sample sizes when T -> Inf
  # d converges to the weighted mean 
  # over the probes, weighted by probe variances	 
  if (nrow(St) == 1) {
    v <- St[1,]
  } else {
    v <- colSums(St / s2) / sum(1 / s2)
  }

  v

}

