

# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


RPA.preprocess <- function (abatch, 
                            bg.method = "rma",
                            normalization.method = "quantiles.robust",
                            cdf = NULL, cel.files = NULL, cel.path = NULL)
{

  # Getting affybatch

  if (is.null(abatch) && (!is.null(cel.files) || !is.null(cel.path))) {
    if (is.null(cel.files) && !is.null(cel.path)) {
       cel.files <- list.celfiles(cel.path, full.names = TRUE)
    }
    abatch <- ReadAffy(filenames = cel.files, compress=getOption("BioC")$affy$compress.cel)  
  } else if (is.null(abatch)) {
    stop("Provide abatch, cel.files or cel.path!")
  } 

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
  abatch <- normalize(abatch, method = normalization.method)
  
  # Log transformation
  message("Logging PM values...")
  q <- log2(pm(abatch))

  message("Retrieving probe positions...")
  # The indices correspond to the output of pm()
  pN <- probeNames(abatch)
  set.inds <- split(1:length(pN), pN) # pNList
  
  return(list(q = q, set.inds = set.inds, cdf = cdf))

}

