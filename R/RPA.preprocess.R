

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
                            cdf = NULL)
{

  message("Preprocessing affybatch...")

  # Set alternative CDF environment if given
  if (!is.null(cdf)) {
	abatch@cdfName <- cdf
	message(paste("Setting alternative CDF", cdf))
  }

  message("Background correcting...")
  abatch2 <- bg.correct(abatch, bg.method, destructive = TRUE)

  message("Normalizing...") 
  abatch2 <- normalize(abatch2, method = normalization.method)
  
  # Log transformation
  message("Logging PM values...")
  q <- log2(pm(abatch2))

  message("Retrieving probe positions...")
  # The indices correspond to the output of pm()
  pN <- probeNames(abatch2)
  set.inds <- split(1:length(pN), pN) # pNList
  
  return(list(q = q, set.inds = set.inds, cdf = cdf))

}
