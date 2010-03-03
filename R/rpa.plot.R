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

rpa.plot <- function (set, rpa.object, highlight.probes = NULL, pcol = "darkgrey", dcol = "black", cex.lab = 1.5, cex.axis = 1) {

  # get the associated affybatch
  abatch <- rpa.object$abatch
  
  # Use alternative CDF environment if given                                        
  if (!is.null(rpa.object$cdf)) { abatch@cdfName <- rpa.object$cdf }

     d <- rpa.object[[set]]$d
    s2 <- rpa.object[[set]]$sigma2
  cind <- rpa.object[[set]]$cind
   dat <- rpa.object$data
  
  # Find probe (pm) indices for this set
  #pmindices <- pmindex(abatch, set)[[1]]
  pN <- probeNames(abatch, set)
  set.inds <- split(1:length(pN), pN) # pNList
  pmindices <- set.inds[[set]]
  
  # Get probes x chips matrix of probe-wise fold-changes
  S <- dat[pmindices, ]
  
  # number of probes
  Np <- length(pmindices)

  # image limits
  ylims <- range(c(as.vector(S), d))

  par(mfrow = c(2, 1))

  # expression figure
  plot(c(1,2,3), type = 'n',
       xlim = c(1, ncol(S)),
       ylim = ylims,
       xlab = "Samples",
       ylab = "Signal log-ratio",
       cex.lab = cex.lab,
       cex.axis = cex.axis,
       main = paste(set,"/ Probe-level values and estimated d"),
       las = 1, xaxt = 'n')
  axis(1, at = 1:3, labels = colnames(rpa.object$data))

  for (i in 1:Np) {
    lines(S[i, ], col = pcol, lwd = 2)
  }

  if (!is.null(highlight.probes)) {
    lines(S[highlight.probes, ], lty = 2, lwd = 2) 
  }

  lines(d, col = dcol, lwd = 2)

  # probe variances
  barplot(s2, main = paste(set,"/ Probe-specific variances (sigma2)"),
          ylab = "Variance", xlab = "Probe index",
          #names.arg = paste(set,seq(length(s2))), las = 2)
          names.arg = 1:length(s2),
          las = 1,
          cex.lab = cex.lab,
          cex.axis = cex.axis)

}
