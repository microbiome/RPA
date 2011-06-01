# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


update.hyperparameters <- function (R, alpha, beta) {
  
  # R: arrays x probes matrix; observations vs. estimated real signal
  #    with mean (d) removed: R <- S - d.

  # return updated hyperparameters
  list(alpha = update.alpha(nrow(R), alpha), beta = update.beta(R, beta))
}


update.alpha <- function (T, alpha) { alpha + T/2 }


update.beta <- function (R, beta, mode = "robust") {

  # FIXME: define separate funcs for modes to speed up?

  if (mode == "approx") {
    # ML estimate, see for instance Bishop p. 100, eqs. 2.150-2.151
    # of variance ((sum(x-x.mean)^2)/n)
    # Note R already assumed to be zero-mean (S = d + N(0, s2))
    # This is same as apply(R, 2, var) except we use 1/n while var uses 1/(n-1)
    # Note this approximation ignores marginalization over the reference sample
    beta + colSums(R^2)/2 
  } else if (mode == "robust") {
    # update beta with posterior mean
    # assuming alpha, beta(hat) are given
    beta.c(beta, R)
  }
}


# Provide compiled version of betahat, about 1.5-fold speedup seen
betahat.f <- function (beta, R) {
     r.colsums <- colSums(R)
    r2.colsums <- colSums(R^2)
    beta + 0.5*(r2.colsums - r.colsums^2 / (nrow(R) + 1))
    # NOTE: by definition,  r.colsums ~ 0 by expectation
    # since R = S - d etc. Therefore we have approximation
    # beta + r2.colsums/2 
    # beta + (N/2)*(r2.colsums/N) 
    # beta + (N/2)* variance(R)
    # vrt. (N/2)*var = 0.5*(N * var) = 0.5 * sum(s^2)
}

# NOTE: this essentially gives 
# beta + (nrow(S)/2) * (ML estimate of (probe)variance);
require(compiler)
beta.c <- cmpfun(betahat.f) 



##################################################################

#(r2.colsums - r.colsums^2 / ((nrow(R) + 1)))  *(1/nrow(R))
#r.colsums <- colSums(R)
#r2.colsums <- colSums(R^2)
#beta + (r2.colsums - r.colsums^2 / ((nrow(R) + 1)))/2
#N <- nrow(R)
#print("1-")
#print((N/2)*colSums(R^2)/(N-1))
#print("2-")
#print((N/2)*apply(R,2,var))
#print("3-")
#print(colSums(R^2)/2)
#print("3-")






############################






