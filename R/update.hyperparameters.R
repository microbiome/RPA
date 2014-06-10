# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2010-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' hyperparameter.update
#' Update hyperparameters Update shape (alpha) and scale (beta) parameters of the inverse gamma distribution.
#'
#' @param dat A probes x samples matrix (probeset).
#' @param alpha Shape parameter of inverse gamma density for the probe variances.
#' @param beta Scale parameter of inverse gamma density for the probe variances.
#' @param th Convergence threshold.
#'
#' @details Shape update: alpha <- alpha + T/2; Scale update: beta <- alpha * s2 where s2 is the updated variance for each probe (the mode of variances is given by beta/alpha). The variances (s2) are updated by EM type algorithm, see s2.update.
#'
#'@return A list with elements alpha, beta (corresponding to the shape and scale parameters of inverse gamma distribution, respectively).
#'
#' @seealso s2.update, rpa.online
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#'## Generate and fit toydata, learn hyperparameters
#'#set.seed(11122)
#'#P <- 11   # number of probes
#'#N <- 5000 # number of arrays
#'#real <- sample.probeset(P = P, n = N, shape = 3, scale = 1, mu.real = 4)
#'#dat <- real$dat # probes x samples#
#'#
#'## Set priors
#'#alpha <- 1e-2
#'#beta  <- rep(1e-2, P)
#'## Operate in batches
#'#step <- 1000
#'#for (ni in seq(1, N, step)) {
#'#  batch <- ni:(ni+step-1)  
#'#  hp <- hyperparameter.update(dat[,batch], alpha, beta, th = 1e-2)
#'#  alpha <- hp$alpha
#'#  beta <- hp$beta
#'#}
#'## Final variance estimate
#'#s2 <- beta/alpha
#'#
#'## Compare real and estimated variances
#'#plot(sqrt(real$tau2), sqrt(s2), main = cor(sqrt(real$tau2), sqrt(s2))); abline(0,1)
#'
#' @export
#' @keywords utilities

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

  # FIXME: mode = "approx" into default if considerable speedups without remarkable compromises

  if (mode == "approx") {
    # ML estimate, see for instance Bishop p. 100, eqs. 2.150-2.151
    # of variance ((sum(x-x.mean)^2)/n)
    # Note R already assumed to be zero-mean (S = d + N(0, s2))
    # This is same as apply(R, 2, var) except we use 1/n while var uses 1/(n-1)
    # Note: this approximation ignores marginalization over the reference sample
    beta + colSums(R^2)/2 
  } else if (mode == "robust") {
    # update beta with posterior mean
    # assuming alpha, beta(hat) are given
    betahat.f(beta, R)
  }
}

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

s2.update <- function (dat, alpha = 1e-2, beta = 1e-2, s2.init = NULL, th = 1e-2, maxloop = 1e6) {

  if (is.null(s2.init)) { s2.init <- rep(1, length(beta)) }

  s2 <- s2.init 
  epsilon <- Inf
  n <- ncol(dat)
  ndot <- n-1
  #n.sqrt <- sqrt(n)
  
  loopcnt <- 0

  while (epsilon > th && loopcnt < maxloop) {

    s2.old <- s2

    # Generalized EM; should be sufficient to improve by using the
    # mode of d (mu.hat) instead of sampling from d (d is a long
    # vector anyway, giving sufficiently many instances); updating with
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

    loopcnt <- loopcnt + 1

  }

  s2

}

s2hat <- function (s2) { 1 / sum(1 / s2) }


s2.neglogp <- function (s2, ndot, alpha, beta.inv, k) {

  # Remove sign; s2 always >0
  s2 <- abs(s2)			  

  # mode of s2 ((N-1) * s2 / s2.obs ~ chisq_(N-1) and mode is df - 2, here df = N-1)
  log.likelihood <- dchisq(s2 * k, df = ndot, log = TRUE)
  log.prior <- dgamma(1/s2, shape = alpha, scale = beta.inv, log = TRUE) # beta.inv <- 1/beta

  # minimize -logP; separately for each probe likelihood + prior
  -(sum(log.likelihood + log.prior))

}


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

# Provide compiled version of approximate beta update
beta.fast <- function (beta, R) {
  beta + colSums(R^2)/2
}

