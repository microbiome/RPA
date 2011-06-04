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


  # The two modes give essentially same results when 
  # sample size is sufficient
  # R <- matrix(rnorm(1000), 100,10)
  # r.colsums <- colSums(R)
  # r2.colsums <- colSums(R^2)
  # plot(colSums(R^2)/2, 0.5*(r2.colsums - r.colsums^2 / (nrow(R) + 1)))
  # abline(0,1)
    # NOTE: by definition,  r.colsums ~ 0 by expectation
    # since R = S - d etc. Therefore we have approximation
    # beta + r2.colsums/2 
    # beta + (N/2)*(r2.colsums/N) 
    # beta + (N/2)* variance(R)
    # vrt. (N/2)*var = 0.5*(N * var) = 0.5 * sum(s^2)

#######################################

# EM update for beta
# P(s2) ~ Prod_s (integrate_d P(s2_p | d) P(d | sigmahat))
# where sigmahat is estimated from previous iteration
# through P(d|sigmahat) ~ Prod_p (P(d_s | x_ps, s2_old_p))
# and P(s2_p|d) ~ Prod_s P(d_s|x_ps, s2_p)

beta.EM <- function (dat, s2.old) {
  
# TODO FIXME TODOTODO


  # dat: probes x samples

  # Prior for d from previous step:
  # P(d | muhat, sigmahat) = Prod_s P(d_s | muhat, sigmahat)
  #  ~ Prod_s Prod_p N(d_s | x_ps, s2__old_p) 
  #  ~ Prod_s N(d_s | muhat_s, s2hat)
  # from product of Gaussians we get 
  # mean: weighted average over the x's
  muhat <- d.update.fast(dat, s2.old)
  # variance (speedup by combining with d.update.fast
  # which also calculates the denominator
  s2hat <- 1/sum(1/s2.old)
  # This is for d.s ~ N(d.s | muhat, s2hat)

  # Now combine with update part (s2 to be estimated):
  # update each probe p separately: P(d.s | x.ps, s2.p)
  # again, we get from product of normal densities
  # P(s2.p) ~ P(d.s | x.ps, s2.p) P(d.s | muhat.s, s2hat)
  muhat.B <- NULL
  for (probe in 1:nrow(dat)) {
    mub <- d.update.fast(rbind(dat[probe,], muhat), c(s2.old[[probe]], s2hat[[probe]]))
    muhat.B <- rbind(muhat.B, mub)
  }
  s2hat.B <- 1/(1/s2.old + 1/s2hat)


}

########################################

# Provide compiled version of approximate beta update
betahat.appr <- function (beta, R) {
  beta + colSums(R^2)/2
}
beta.fast.c <- cmpfun(betahat.appr) 

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






