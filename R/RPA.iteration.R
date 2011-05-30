
# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


# Changelog:

RPA.iteration <- function(S,
                          epsilon = 1e-3,
                            alpha = NULL,
                             beta = NULL,
                    sigma2.method = "robust",
                         d.method = "fast",
                          maxloop = 1e6)
{


  P <- ncol(S) # number of probes
  T <- nrow(S) # Number of arrays (except reference)

  # Check: if affybatch/probeset is erroneous and contains just NAs or NaNs then return NA vector
  if (all(is.nan(S) | is.na(S))) { 
    return(list(d = rep(NA, T), sigma2 = rep(NA, P))) 
  }

  # not used in computation,
  # just to check convergence at first iteration
  d <- rowMeans(S) # initialize by mean over the probes
  d.old <- (-d) 
  
  # uninformative priors for sigma2.methods mean, mode, var;
  # informative for 'robust', or alpha, beta are provided by user
  alpha <- set.alpha(alpha, sigma2.method, P)
  beta <- set.beta(alpha, sigma2.method, P)

  # Update hyperparameters
  alpha <- update.alpha(T, alpha) # FIXME: if online version turns out problematic, check this:
  	   		   	  # is it ok to update alpha already here?
  # Do NOT update beta yet; it will be updated in while loop
  # after estimating d based on the current priors!

  # Confirm that sigma2.method is valid for these parameters
  if (sigma2.method == "mean" || sigma2.method == "robust") {
    ifelse(all(alpha > 1), TRUE, stop("alpha > 1 - (N.arrays - 1) / 2 required for this sigma2.method"))
  } else {}

  # initialize sigma2
  sigma2 <- RPA.sigma2.update(R, alpha, beta, sigma2.method)
  #sigma2 <- rep.int(1, P)

  # optimize until convergence
  loopcnt <- 0

  if (d.method == "fast") {
    while ((max(abs(c(d - d.old))) > epsilon) && loopcnt < maxloop) {

      d.old <- d

      # update d, given sigma2
      d <- d.update.fast(t(S), sigma2)

      # Estimate noise 
      R <- S - d

      # beta update
      beta <- update.beta(R, beta)

      # update sigma2
      sigma2 <- RPA.sigma2.update(R, alpha, beta, sigma2.method)

      # follow iteration count to avoid potentially infinite loops
      loopcnt <- loopcnt + 1 
    }
  } else if (d.method == "basic") {

      while ((max(abs(c(d - d.old))) > epsilon) && loopcnt < maxloop) {

        # separate while loops for d.methods to avoid logical comparisons 
        # during iteration	
        d.old <- d

        # beta update
        beta <- update.beta(S - d, beta)

        # update d
        d <- optim(d, fn = RPA.dcost, method = "BFGS", sigma2 = sigma2, S = S)$par

        # update sigma2
        sigma2 <- RPA.sigma2.update(S - d, alpha, beta, sigma2.method)

        # follow iteration count to avoid potentially infinite loops
        loopcnt <- loopcnt + 1 

      }
   } 

  list(d = d, sigma2 = sigma2, alpha = alpha, beta = beta)
}
