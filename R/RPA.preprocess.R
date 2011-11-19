

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

get.set.inds <- function (cel.files, cdf = NULL, sets = NULL) {

  # Get probe position indices
  abatch <- ReadAffy(filenames = cel.files[1:2], compress=getOption("BioC")$affy$compress.cel)
  # Set alternative CDF environment if given
  if (!is.null(cdf)) { abatch@cdfName <- cdf }    

  # Check names for the investigated probesets
  # if my.sets not specified, take all sets in abatch
  if ( is.null(sets) ) { sets <- geneNames(abatch) } 

  # Retrieve probe positions
  # The indices correspond to the output of pm()
  pN <- probeNames(abatch)
  set.inds <- split(1:length(pN), pN)[sets] # pNList
  names(set.inds) <- sets
  
  #list(set.inds = set.inds, cdf = cdf)
  set.inds
 
}

get.batches <- function (items, batch.size, shuffle = FALSE) {

  
  # Random ordering for the items?
  if (shuffle) {items <- sample(items)}

  # N elements into batches of size batch.size
  # last batch can be smaller

  if (length(items) == 1 && is.numeric(items)) {
    N <- items
    items <- 1:N
  } else {
    N <- length(items)
  }

  if (N < batch.size) {
     warning("batch.size > N, setting batch.size = N.")
  }

  batches <- list()
  ns <- floor(c(0, seq(batch.size, N, batch.size)))
  if (ns[[length(ns)]] < N) {ns[[length(ns) + 1]] <- N}
  cnt <- 0
  for (i in 2:length(ns)) {
    cnt <- cnt + 1
    n.start <- ns[[i-1]]
    n.stop <- ns[[i]]
    batches[[i-1]] <- items[(n.start+1):n.stop]
  }

  # Provide ID for each batch
  names(batches) <- paste("batch", 1:length(batches), sep = "-")
  
  batches
}

