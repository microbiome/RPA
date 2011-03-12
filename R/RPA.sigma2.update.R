# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


RPA.sigma2.update <- function (d, S = S, alphahat, beta, sigma2.method = "robust") {

  #
  # S: arrays x probes matrix; observations vs. estimated real signal
  #
  
  R <- S - d 

  if (!sigma2.method == "var") {

    # update betahat with posterior mean
    # assuming alpha, beta(hat) are given
    r.colsums <- colSums(R)
    r2.colsums <- colSums(R^2)
    betahat <- .5*(2 * beta + r2.colsums - r.colsums^2 / (nrow(S) + 1))

    # update betahat by posterior mean (sigma2.method 'mean' and
    # 'robust') or mode ('mode')
    # alphahat > 1 required in sigma2.method 'mode' (and 'robust')
    # betahat / (alphahat - 1) = sigma2 mean for invgam(sig^2 | alphahat,betahat)
    # betahat / (alphahat + 1) = sigma2 mode for invgam(sig^2 | alphahat,betahat)
    # (sigma2.method == "mean" || sigma2.method == "robust")
    if (!sigma2.method == "mode") {
      s2 <- betahat / (alphahat - 1)
    } else {s2 <- betahat / (alphahat + 1)}
  } else if (sigma2.method == "var") {
    # Assume uninformative priors alpha, beta -> 0	  
    # NOTE: RPA converges to variance with large sample size
    s2 <- apply(R, 2, var) 		
  } else { stop("sigma2.method missing!") }

  # return updated sigma2
  s2
}

