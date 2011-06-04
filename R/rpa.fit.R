# This file is a part of the RPA (Robust Probabilistic Averaging)
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 
rpa.fit <- function (dat, cind = 1, epsilon = 1e-2, alpha = NULL, beta = NULL, sigma2.method = "fast", d.method = "fast") {

  # dat: original data (probes x samples)
  
  # Fits RPA on fold-change data calculated against the reference sample.
  # After estimating the RPA fold-changes, also fits the mean in the
  # original data domain since this is often desired.

  # Extract reference sample  
  # Get samples x probes matrix of probe-wise fold-changes
  if (is.null(colnames(dat))) {colnames(dat) <- 1:ncol(dat)}

  # Fit RPA
  S <- t(dat[, -cind] - dat[, cind])
  estimated <- RPA.iteration(S, epsilon, alpha, beta, sigma2.method, d.method)

  # Estimate overall signal
  if (d.method == "fast") {
    mu <- d.update.fast(dat, estimated$sigma2)
    affinity <- estimate.affinities(dat, mu)
    mu.abs <- mean(mu[-cind] - estimated$d) # Difference between total signal and shape variable
  } else if (d.method == "basic") {
    # add the reference
    d <- rep.int(0, ncol(dat))
    names(d) <- colnames(dat)
    d[rownames(S)] <- estimated$d
    prs <- t(dat) - d
    # overall signal level is obtained by weighted mean
    # weighted by probe-specific affinities
    # finally take mean over all arrays to get robust estimate
    # (all rows should give approximately same result)
    mu.abs <- mean(rowSums(prs/estimated$sigma2)/sum(1/estimated$sigma2))
    # prs = mu.abs + af + noise; dat = d + mu.abs + af + noise
    affinity <- colMeans(prs) - mu.abs
    mu <- mu.abs + d

  }

  # Now return final fitted parameters
  # d: differential signal between reference sample and other samples
  # mu: fitted signal in original data domain
  # affinity: probe-specific affinities
  # sigma2: probe-specific stochastic noise
  # model for probe p: x_p = mu.orig + affinity_p + noise_p; noise ~ N(0, sigma2_p)

  new("rpa.fit",
      list(mu = mu, mu.real = mu.abs, 
           sigma2 = estimated$sigma2, 
           affinity = affinity, 
           data = dat,
	   alpha = estimated$alpha,
           beta = estimated$beta))
  
}
