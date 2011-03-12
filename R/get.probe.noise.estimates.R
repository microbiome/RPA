
# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


get.probe.noise.estimates <- function (rpa.res, sets = NULL, normalization = NULL, verbose = FALSE) {

  if (is.null(sets)) {sets <- names(rpa.res)}
  
  if ( is.null( normalization ) ) {
    # No normalization: return probe variances
    rpa.res$sigma2 <- rpa.res$sigma2

  } else if (normalization == "withinset.categorical") {

    for (set in sets) {
      # Order probes within the probeset by their noise level
      # smallest index means least noise
      if (verbose) {cat(set)}
      o <- order(rpa.res$sigma2[[set]])
      rpa.res$sigma2[[set]][o] <- 1:length(o)      
    }
    
  } else if (normalization == "withinset.relative") {

    for (set in sets) {
      # compare probe standard deviation to probeset signal standard
      # deviation
      # within the probeset
      # smallest means least noise
      if (verbose) {cat(set)}
      vars <- rpa.res$sigma2[[set]]
      rpa.res$sigma2[[set]] <- sqrt(vars)/sqrt(var(rpa.res$d[set,]))
    }
  } else if (normalization == "withinset.weights") {

    for (set in sets) {
      # inverse of the weights used to compute
      # probeset-level summary signal d
      # within the probeset
      # smallest means least noise (largest weight = smallest inverse weight)
      if (verbose) {cat(set)}
      vars <- rpa.res$sigma2[[set]]
      rpa.res$sigma2[[set]] <- 1/((1 / vars) / sum(1 / vars))
    }
  }

  rpa.res$sigma2[sets]

}
