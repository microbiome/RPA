# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

estimate.hyperparameters <- function (priors, set.inds, batches, cdf, quantile.basis, bg.method, normalization.method, epsilon, cind, load.batches = NULL, mc.cores = 1, verbose = TRUE) {

  #load.batches = batch.file.id; mc.cores = 4
  
  # Hyperparameter estimation through batches

  # Initialize hyperparameters

  # Note: alpha is scalar and same for all probesets 
  # alpha <- alpha + N/2 at each batch

  sets <- names(set.inds)
  
  if (verbose) {message("Initialize priors")}
  alpha <- priors$alpha # initialize 
  betas <- mclapply(sets, function (set) { rep(priors$beta, length(set.inds[[set]])) }, mc.cores = mc.cores)
  names(betas) <- sets

  for (i in 1:length(batches)) {

    if (verbose) {message(paste("Updating hyperparameters on batch", i, "/", length(batches)))}
    
    # Get background corrected, quantile normalized, and logged probe-level matrix
    if (!is.null(load.batches)) {
      batch.file <- paste(load.batches, "-", i, ".RData", sep = "")
    }
      
    q <- get.probe.matrix(cels = batches[[i]], cdf, quantile.basis, bg.method, normalization.method, batch.file, verbose = verbose)

    # Get probes x samples matrix of probe-wise fold-changes
    # q <- matrix(q[, -cind] - q[, cind], nrow(q))
    T <- ncol(q) # Number of arrays expect reference

    # Get probes x samples matrices for each probeset

    q <- mclapply(set.inds, function (pmis) { matrix(q[pmis,], length(pmis)) }, mc.cores = mc.cores)
    names(q) <- sets	    

    # Update variance for each probeset
    s2s <- mclapply(sets, function (set) {
      s2.update(q[[set]], alpha, betas[[set]], s2.init = betas[[set]]/alpha, th = epsilon)
    }, mc.cores = mc.cores) # FIXME move conv. param. to arguments
    names(s2s) <- sets
    
    # Update alpha, beta (variance = beta/alpha at mode with large T)
    alpha <- update.alpha(T, alpha)
    betas <- mclapply(s2s, function (s2) { s2 * alpha }, mc.cores = mc.cores)

  }

  # Get final estimated variances for each probeset based on hyperparameter posteriors
  variances <- mclapply(betas, function (beta) {beta/alpha}, mc.cores = mc.cores)
  names(variances) <- names(betas) 

  
  list(alpha = alpha, betas = betas, variances = variances)  

}

hyperparameter.update <- function (dat, alpha, beta, th = 1e-2) {
  
  # dat: probes x samples matrix

  # Estimate s2 with EM-type procedure
  s2 <- s2.update(dat, alpha, beta, s2.init = beta/alpha, th = th)

  # Update hyperparams: 
  alpha <- update.alpha(ncol(dat), alpha)
  #beta <- alpha * s2 # from the mode of s2, i.e. beta/alpha

  # return updated hyperparameters
  list(alpha = alpha, beta = alpha * s2)
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
    betahat.f(beta, R)
  }
}


s2.update <- function (dat, alpha = 1e-2, beta = 1e-2, s2.init = NULL, th = 1e-2) {

  if (is.null(s2.init)) { s2.init <- rep(1, length(beta)) }

  s2 <- s2.init 
  epsilon <- Inf
  n <- ncol(dat)
  ndot <- n-1
  #n.sqrt <- sqrt(n)
  
  while (epsilon > th) {

    s2.old <- s2

    # Generalized EM; should be sufficient to improve by using the
    # mode of d (mu.hat) instead of sampling from d (d is a long
    # vector anyway, giving sufficiently many instants); updating with
    # the mode here as complete sampling would slow down computation
    # considerably, probably without much gain in performance in this
    # case.
    d <- d.update.fast(dat, s2)
    s2hat <- s2hat(s2)
    
    # optimize s2; regularize observed variance by adding small
    # constant on observed variances; use here s2hat as an adaptive
    # lower bound (s2hat is variance of the mean). The regularization
    # is needed as with the current implementation (which uses mode of
    # d instead of full sampling over d space) optimization gets
    # sometimes stuck to local optima, collapsing s2 -> 0 for one of
    # the probes (about 1e-4~1e-5 occurrence frequency); the prior
    # does not seem to prevent this here.
    s2.obs <- s2obs(dat, d, ndot) + s2hat
     k <- ndot / s2.obs
    s2 <- abs(optim(s2, fn = s2.neglogp, ndot = ndot, alpha = alpha, beta.inv = 1/beta, k = k, method = "BFGS")$par)

    # FIXME: check if 'optimize' would be faster?
    epsilon <- max(abs(s2 - s2.old))

  }

  s2

}

# FIXME: could be utilized in d.update.fast.c to speed up
s2hat <- function (s2) { 1 / sum(1 / s2) }
#s2hat.c <- cmpfun(s2hat) 

s2.neglogp <- function (s2, ndot, alpha, beta.inv, k) {

  # Remove sign; s2 always >0
  s2 <- abs(s2)			  

  # mode of s2 ((N-1) * s2 / s2.obs ~ chisq_(N-1) and mode is df - 2, here df = N-1)
  log.likelihood <- dchisq(s2 * k, df = ndot, log = TRUE)
  log.prior <- dgamma(1/s2, shape = alpha, scale = beta.inv, log = TRUE) # beta.inv <- 1/beta

  # minimize -logP; separately for each probe likelihood + prior
  -(sum(log.likelihood + log.prior))

}


#dchisq.c <- cmpfun(dchisq) 
#dgamma.c <- cmpfun(dgamma) 

s2obs <- function (dat, d, ndot) {
  # Center the probes to obtain approximation:
  # x = d + mu.abs + mu.affinity + epsilon.reference + epsilon.samples
  # by centering we remove mu.abs + mu.affinity + epsilon.reference
  # with large sample size it is safe to assume that epsilon.reference = 0
  # so this does not affect the results, compared to the exact solution
  # where the variance in epsilon is estimated also including epsilon.reference ie. 
  # x = d + epsilon.reference + epsilon.samples
  # which would just shift x little but as it is marginalized in the treatment the
  # result will converge to estimating variance from 
  # x = d + epsilon.samples when N -> Inf. Therefore, we use:
  # datc <- t(centerData(t(dat)))

  datc <- centerData(t(dat) - d)

  # avoid getting stuck to local optima with d (occurs sometimes ~
  # 1e-4) by adding small constant on observed variance (which should
  # never be 0 in real data) seems that the prior term does not always
  # prevent such collapse with the current implementation; can
  # probably be improved so that regularization is cpompletely handled
  # by the prior (FIXME)
  #s2 + 1e-3
  colSums(datc^2)/ndot
}

#s2obs.c <- cmpfun(s2obs) 


########################################

# Provide compiled version of approximate beta update
beta.fast <- function (beta, R) {
  beta + colSums(R^2)/2
}
#beta.fast.c <- cmpfun(betahat.fast) 

#######################################

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
#require(compiler)

#beta.c <- cmpfun(betahat.f) 










