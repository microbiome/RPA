

# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


RPA.preprocess <- function (abatch, cind = 1,
                            bg.method = "rma",
                            normalization.method = "quantiles.robust",
                            cdf = NULL, quantile.n = 50)
{

  message("Preprocessing affybatch...")

  # Set alternative CDF environment if given
  if (!is.null(cdf)) {
	abatch@cdfName <- cdf
	message(paste("Setting alternative CDF", cdf))
  }

  message("Background correcting...")
  abatch <- bg.correct(abatch, bg.method, destructive = TRUE)
  # FIXME: here the abatch values for some reason are set to NaNs!
  # for defected affybatch
  
  message("Normalizing...") 
  if (normalization.method == "quantiles.online") {
    abatch <- online.quantile(abatch, quantile.n)
  } else {
    abatch <- normalize(abatch, method = normalization.method)
  }

  # Log transformation
  message("Logging PM values...")
  q <- log2(pm(abatch))

  message("Retrieving probe positions...")
  # The indices correspond to the output of pm()
  pN <- probeNames(abatch)
  set.inds <- split(1:length(pN), pN) # pNList
  
  return(list(q = q, set.inds = set.inds, cdf = cdf))

}
