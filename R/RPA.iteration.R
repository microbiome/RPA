# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

RPA.iteration <- function(S,
                          epsilon = 1e-3,
                            alpha = NULL,
                             beta = NULL,
                    tau2.method = "fast",
                         d.method = "fast",
                          maxloop = 1e6)
{


  P <- ncol(S) # number of probes
  T <- nrow(S) # Number of arrays (except reference)

  # Accommodate single-probe probesets
  if (P == 1) { return(list(d = as.vector(S), tau2 = 0, alpha = T/2, beta = 0)) }

  # Check: if affybatch/probeset is erroneous and contains just NAs or NaNs then return NA vector
  if (all(is.nan(S) | is.na(S))) { 
    return(list(d = rep(NA, T), tau2 = rep(NA, P))) 
  }

  # uninformative priors for tau2.methods mean, mode, var;
  # informative for 'robust', or alpha, beta are provided by user
  alpha.prior <- alpha <- set.alpha(alpha, tau2.method, P)
  beta.prior <- beta <- set.beta(beta, tau2.method, P)

  # Confirm that alpha is valid for tau2.method 
  if (tau2.method == "mean" || tau2.method == "robust") {
    ifelse(all(alpha > 1), TRUE, stop("alpha > 1 - (N.arrays - 1) / 2 required for this tau2.method"))
  } else {}

  ###############################

  # initialize tau2 with user-defined priors

  if (tau2.method == "var") {
    s2.meth <- "mean"
  } else {
    s2.meth <- tau2.method
  }

  # FIXME: add chance to give user-defined priors in here!
  tau2 <- RPA.tau2.update(NULL, alpha.prior, beta.prior, s2.meth)
  # initialize with large initial variance if not 
  if (length(tau2) == 0) {tau2 <- rep(3*sd(as.vector(S)), P)} 

  # check convergence at first iteration, 
  # not used in calculations
  tau2.old <- rep(Inf, P)
  
  ###############################

  # Update alpha
  # Do NOT update beta yet; it will be updated in while loop
  # after estimating d based on the current priors!
  alpha <- update.alpha(T, alpha) 

  #################################

  # optimize until convergence
  loopcnt <- 0

  if (d.method == "fast") {

    while ((max(abs(c(tau2 - tau2.old))) > epsilon) && loopcnt < maxloop) {

      tau2.old <- tau2

      # update d, given tau2
      d <- d.update.fast(t(S), tau2)

      # Estimate noise 
      R <- S - d

      # beta update (feed in beta prior, not updates from this loop!)
      beta <- update.beta(R, beta.prior)

      # update tau2
      tau2 <- RPA.tau2.update(R, alpha, beta, tau2.method)

      # follow iteration count to avoid potentially infinite loops
      loopcnt <- loopcnt + 1 

    }

  } else if (d.method == "basic") {

      # initialize d for the first iteration
      d <- rowMeans(S)  

      while ((max(abs(c(tau2 - tau2.old))) > epsilon) && loopcnt < maxloop) {

        # separate while loops for d.methods to avoid logical comparisons 
        # during iteration	
	tau2.old <- tau2

        # update d
        d <- optim(d, fn = RPA.dcost, method = "BFGS", tau2 = tau2, S = S)$par

        # Estimate noise 
        R <- S - d

        # beta update (feed in beta prior, not updates from this loop!)
        beta <- update.beta(R, beta.prior)

        # update tau2
        tau2 <- RPA.tau2.update(R, alpha, beta, tau2.method)

        # follow iteration count to avoid potentially infinite loops
        loopcnt <- loopcnt + 1 

      }
   } 

  list(d = d, tau2 = tau2, alpha = alpha, beta = beta)
}

#################################################################################

# Changelog:

RPA.iteration.fast <- function(S,
                          epsilon = 1e-3,
                            alpha.prior,
			    alpha.posterior,
                             beta.prior,
                          maxloop = 1e6)
{

  # S: samples x probes 

  # Check: if affybatch/probeset is erroneous and contains just NAs or NaNs then return NA vector
  if (all(is.nan(S) | is.na(S))) { 
    return(list(d = rep(NA, nrow(S)), tau2 = rep(NA, ncol(S)))) 
  }

  # Initialize variances by the original user-specified priors 
  tau2 <- beta.prior/alpha.prior # approx. mean and mode when sample size is largish

  # check convergence at first iteration, 
  # not used in calculations
  tau2.old <- rep(-1e6, length(tau2))
  
  ###############################

  # optimize until convergence
  loopcnt <- 0

  while ((max(abs(c(tau2 - tau2.old))) > epsilon) && loopcnt < maxloop) {

      tau2.old <- tau2

      # update d, given tau2
      d <- d.update.fast(t(S), tau2)

      # Estimate noise 
      R <- S - d

      # beta update (feed in original beta prior, not updates from this loop!)
      beta.posterior <- beta.fast(beta.prior, R) 
      #betahat.f(beta.prior, R) # update.beta(R, beta.prior)

      # update tau2 / NOTE: use posterior alpha and beta here!
      tau2 <- beta.posterior/alpha.posterior
      #RPA.tau2.update(R, alpha.posterior, beta.posterior, tau2.method = fast)

      # follow iteration count to avoid potentially infinite loops
      loopcnt <- loopcnt + 1 

  }


  list(d = d, tau2 = tau2, beta = beta.posterior)
}
