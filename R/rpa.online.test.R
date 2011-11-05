# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  
rpa.online.test <- function (
		dat,
                  priors = list(alpha = 2, beta = 1),
                 epsilon = 1e-2, 
                    cind = 1,
                 verbose = FALSE,
               bg.method = "rma",
                     cdf = NULL, 
              batch.size = 10) 
{

  batches <- get.batches(1:ncol(dat), batch.size)

  # Initialize hyperparameters
  alpha <- priors$alpha # initialize 
  beta <- rep(priors$beta, nrow(dat))

  for (i in 1:length(batches)) {

    my.dat <- dat[, batches[[i]]]

    # Get probes x samples matrix of probe-wise fold-changes
    #q <- matrix(q[, -cind] - q[, cind], nrow(q))
    T <- ncol(my.dat) # Number of arrays expect reference

    # Update hyperparams: 
    alpha.old <- alpha
    alpha <- update.alpha(ncol(my.dat), alpha)
 
    s2 <- s2.update(my.dat, alpha.old, beta, s2.init = beta/alpha, th = 1e-2)

    beta <- alpha * s2 # from the mode of s2, i.e. beta/alpha

  }

  # Now hyperparameters for probe variances have been estimated
  # Get probeset-level signal estimate by online sweep with the fixed variances

  dest <- c()
  for (i in 1:length(batches)) {

    q <- dat[, batches[[i]]]

    # Get estimated variances for each probeset based on hyperparameter posteriors
    variances <- beta/alpha

    dest <- c(dest, d.update.fast(q, variances))

  }

  list(alpha = alpha, beta = beta, d = dest)
  
}


