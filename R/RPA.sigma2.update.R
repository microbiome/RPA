# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

RPA.sigma2.update <- function (R, alpha, beta, sigma2.method = "robust") {

  # FIXME: online mode not necessarily needed at all here
  # FIXME: speedup by defining the method outside this looping function

  # sigma2 values are notably smaller with 
  # sigma2.method = "fast" than sigma2.method = "robust"
  # This is because different priors for alpha are used
  # they lead to similar probe weights, though so not so much effect on probeset level estimates

  # R <- S - d: arrays x probes matrix; observations vs. estimated real signal
  # Note: alpha here is alphahat = T/2 + alpha w.r.t. user-defined alpha prior

  # alpha > 1 required in sigma2.method 'mode' and 'robust'
  if (sigma2.method == "mean" || sigma2.method == "online" || sigma2.method == "robust") {
    # mean for invgam(sig^2 | alpha,beta)
    s2 <- beta / (alpha - 1) # FIXME: speedup by precalculating alpha - 1?
  } else if (sigma2.method == "mode") {
    # mode for invgam(sig^2 | alpha,beta)
    s2 <- beta / (alpha + 1) 
  } else if (sigma2.method == "fast") {
    # this is often also used
    s2 <- beta / alpha 
  } else if (sigma2.method == "var") {
    # Assume uninformative priors alpha, beta -> 0	  
    # NOTE: RPA converges to variance with large sample size
    # priors not used
    # Do not center: S - d is already assumed to be centered by definition
    # S = d + N(0, sigmaj2)  
    # R <- S - d 
    # faster way to calculate variance than apply(R, 2, var):
    # Note that R is assumed to be centered already,
    # so actually this is more accurate for our case than apply(R,2,var)
    # which would center the data
    #s2 <- colSums(centerData(R)^2)/(nrow(R) - 1) 
    s2 <- colSums(R^2)/(nrow(R) - 1) 
  } else { 
    stop("Invalid sigma2.method provided!") 
  }

  # return updated sigma2
  s2
}


