#
# This file is a part of the RPA program (Robust Probabilistic
# Averaging), see http://www.cis.hut.fi/projects/mi/software/RPA/
#
# Copyright (C) 2008-2010 Leo Lahti (leo.lahti@iki.fi)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License 2 for more details.
# 

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

