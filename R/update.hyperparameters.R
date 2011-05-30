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


update.alpha <- function (T, alpha) { T/2 + alpha }



update.beta <- function (R, beta) {

  # ML estimate, see for instance Bishop p. 100, eqs. 2.150-2.151
  # of variance ((sum(x-x.mean)^2)/n)
  # Note R already assumed to be zero-mean (S = d + N(0, s2))
  # This is same as apply(R, 2, var) except we use 1/n while var uses 1/(n-1)
  beta + colSums(R^2)/2 # shortcut

}


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






