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

RPA.iteration <- function(S,
                          epsilon = 1e-3,
                            alpha = NULL,
                             beta = NULL,
                    sigma2.method = "robust",
                         d.method = "fast",
                          maxloop = 1e6)
{

  # if no prior has been given, use noninformative default priors
  P <- ncol(S) # number of probes
  T <- nrow(S) # Number of arrays (except control)
  S.means <- rowMeans(S)
  St <- t(S)
  
  # initialize with equal weight for all probes
  sigma2 <- rep.int(1, P)

  # uninformative priors for sigma2.methods mean, mode, var;
  # informative for 'robust'

  # if alpha is scalar, set identical prior for all probes with this value
  if (is.null(alpha) && !sigma2.method == "robust") {
    alpha <- rep.int(1e-6, P)
  } else if (is.null(alpha) && sigma2.method == "robust") {
    alpha <- rep.int(2, P)
  } else if (length(alpha) == 1) {
    alpha <- rep.int(alpha, P)
  } else {}

  # if beta is scalar, set identical prior for all probes with this value
  if (is.null(beta) && !sigma2.method == "robust") {
    beta <- rep.int(1e-6, P)
  } else if (is.null(beta) && sigma2.method == "robust") {
    beta <- rep.int(1, P)
  } else if (length(beta) == 1) {
    beta <- rep.int(beta, P)
  } else {}

  # prior (note: set here after all alpha checks)
  alphahat <- T/2 + alpha
  
  # Confirm that sigma2.method is valid for these parameters
  if (sigma2.method == "mean") {
    ifelse(all(alphahat > 1), TRUE, stop("alpha > 1-nrow(S)/2 required for sigma2.method = mean"))
  } else {}
    
  # not used in computation,
  # just to check convergence at first iteration
  sigma2.old <- rep.int(1e3*epsilon, P)

  # These initial values not used in computation,
  # just to check convergence at first iteration
  d <- rep.int(max(S), T) 
  d.old <- (-d) 
  d.init <- rep.int(0,T) 

  # optimize until convergence
  loopcnt <- 0

  if (d.method == "fast") {
    while ((max(abs(c(d - d.old))) > epsilon) && loopcnt < maxloop) {

      d.old <- d
      sigma2.old <- sigma2

      # update d
      d <- d.update.fast(St, sigma2)

      # update sigma2	
      sigma2 <- RPA.sigma2.update(d, S, alphahat, beta, sigma2.method)

      # follow iteration count to avoid potentially infinite loops
      loopcnt <- loopcnt + 1 
    }
   } else if (d.method == "basic") {

      while ((max(abs(c(d - d.old))) > epsilon) && loopcnt < maxloop) {

        # separate while loops for d.methods to avoid logical comparisons 
        # during iteration	
        d.old <- d
        sigma2.old <- sigma2

        # update d
        d <- optim(d.init, fn = RPA.dcost, method = "BFGS", sigma2 = sigma2, S = S)$par

        # update sigma2
        sigma2 <- RPA.sigma2.update(d, S, alphahat, beta, sigma2.method)

        # follow iteration count to avoid potentially infinite loops
        loopcnt <- loopcnt + 1 

      }
   } 

  list(d = d, sigma2 = sigma2)
}
