# This file is a part of the RPA (Robust Probabilistic Averaging)
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 
rpa.fit <- function (dat, cind = 1, epsilon = 1e-2, alpha = NULL, beta = NULL, sigma2.method = "robust", d.method = "fast", affinity.method = "rpa") {

  # dat: original data (probes x samples)
  
  # Fits RPA on fold-change data calculated against the reference sample.
  # After estimating the RPA fold-changes, also fits the mean in the
  # original data domain since this is often desired.

  # Extract reference sample  
  # Get samples x probes matrix of probe-wise fold-changes
  if (is.null(colnames(dat))) {colnames(dat) <- 1:ncol(dat)}
  S <- t(dat[, -cind] - dat[, cind])

  # Fit RPA
  estimated <- RPA.iteration(S, epsilon, alpha, beta, sigma2.method, d.method)

  # Retrieve results. Reference sample is zero by definition.
  # Add it to the result since it was is used in the calculation
  d <- rep.int(0, ncol(dat))
  names(d) <- colnames(dat)
  d[rownames(S)] <- estimated$d

  # Now d gives signal shape.
  # Next: evaluate mean level in the original data
  if (affinity.method == "rpa") {
    # weigh probes by their estimated reliability level also when
    # evaluating affinities. Heuristic solution but expected to be better
    # than giving all probes equal weight in affinity calculation, as in the
    # 'zeromean' option
    sigmas <- estimated$sigma2
  } else if (affinity.method == "zeromean") {
    # require that probe affinities sum to zero
    # (similar criteria as in RMA)
    sigmas <- rep.int(1, length(estimated$sigma2))
  }  

  mu <- estimate.affinities(dat, d, sigmas)

  # Now return final fitted parameters
  # d: differential signal between reference sample and other samples
  # mu: fitted signal in original data domain
  # affinity: probe-specific affinities
  # sigma2: probe-specific stochastic noise
  # model for probe p: x_p = mu.orig + affinity_p + noise_p; noise ~ N(0, sigma2_p)

  new("rpa.fit",
      list(mu = d + mu$real, mu.real = mu$real, 
           sigma2 = estimated$sigma2, 
           affinity = mu$probe, 
           data = dat,
	   alpha = estimated$alpha,
           beta = estimated$beta))
  
}
