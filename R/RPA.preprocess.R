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

RPA.preprocess <- function (abatch, cind = 1,
                            bg.method = "rma",
                            normalization.method = "quantiles.robust",
                            cdf = NULL)
{

  message("Preprocessing affybatch:")

  # Set alternative CDF environment if given
  if (!is.null(cdf)) {
	abatch@cdfName <- cdf
	message(paste("Using alternative CDF environment:", cdf))
  }

  message("Background correcting")
  abatch2 <- bg.correct(abatch, bg.method)

  message("Normalizing") 
  abatch2 <- normalize(abatch2, method = normalization.method)

  # Log transformation
  #pmindex(Dilution, sets[[1]])[[1]]
  #q <- log2(exprs(abatch2)) 
  #fcmat[pmindices, ]
  q <- log2(pm(abatch2))

  # Pick probe name table.
  # The indices correspond to the output of pm()
  pN <- probeNames(abatch2)
  set.inds <- split(1:length(pN), pN) # pNList
  
  message(paste("Setting control array: ", cind, " (", colnames(q)[[cind]], ")",sep=""))
  fcmat <- q[, -cind] - q[, cind]
  colnames(fcmat) <- colnames(q)[-cind]

  return(list(fcmat = fcmat, cind = cind, cdf = cdf, set.inds = set.inds))
}
