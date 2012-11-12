# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  # To estimate means in the original data domain let us assume that
  # each probe-level observation x is of the following form:
  # x = d + a + noise,
  # where x and d are vectors over samples,
  # a is a scalar (vector with identical elements)
  # noise is Gaussian with zero mean and probe-specific variance parameters tau2 
  # Then the parameter a will indicate how much probe-level observation
  # deviates from the estimated signal shape d.
  # This deviation is further decomposed as
  # a = a.real + mu.probe, where
  # a.real describes the 'real' signal level, common for all probes
  # mu.probe describes probe affinity effect
  # Let us now assume that mu.probe ~ N(0, sigma.probe).
  # This encodes the assumption that in general the affinity
  # effect of each probe tends to be close to zero.

  # if identical sigma.probe is used for all probes then a.real is
  # estimated by the average of the probe effects a. Note that then
  # probe-specific affinities mu.probe will sum to exactly zero, which
  # gives an analogous model than used in RMA, which uses this
  # assumption to fit medianpolish.

  # Probes can also be weighted by setting probe-specific sigmas.  In
  # "rpa" option, set sigma.probe to the tau2 value of the probe
  # estimated by RPA.  Note that while tau2 in RPA measures
  # stochastic noise, and NOT the affinity effect, we can use it
  # also to weigh the probes according to how much
  # they contribute to the overall signal shape. Then probes
  # that have little effect on the signal shape (i.e. are very noisy
  # and likely to be contaminated by many unrelated signals) 
  # also contribute less to the absolute signal estimate. If no other
  # prior information is available, using stochastic parameters tau2
  # to determine probe weights is likely to work better than simple
  # averaging of the probes without weights. Also in this case the
  # probe affinities sum close to zero but there is some flexibility,
  # and more noisy probes can be downweighted.

  # Mean level for each probe after extracting the
  # signal shape. This is the ML estimate for mu = mu.real + mu.probe
  # after ignoring the negligible extra variance from the measurement noise
  # averaged over all probes and samples (converges rapidly to zero)
  # a <- colMeans(t(dat) - d)

  # weighted mean, weighted by probe-specific affinity priors
  # a.real <- sum(a/sigmas)/sum(1/sigmas)
  
  # FIXME: add other ways to define affinity priors later

estimate.affinities <- function (dat, a) {

  # In RPA, the final signal estimate
  # is weighted average over the probes, 
  # both in terms of shape and affinities
  
  # Calculate affinities w.r.t. overall signal
  colMeans(t(dat) - a)

}
