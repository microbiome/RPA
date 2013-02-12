# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' RPA.preprocess
#' Preprocess AffyBatch object for RPA.
#'
#' @param abatch An AffyBatch object. 
#' @param bg.method Specify background correction method. See bgcorrect.methods(abatch) for options.
#' @param normalization.method Specify normalization method. See normalize.methods(abatch) for options. For memory-efficient online version, use "quantiles.online".
#' @param cdf The CDF environment used in the analysis. 
#' @param cel.files List of CEL files to preprocess.
#' @param cel.path Path to CEL file directory.
#' @param quantile.basis Optional. Basis for quantile normalization.
#
#' @details Background correction, quantile normalization and log2-transformation for probe-level raw data in abatch. Then probe-level differential expression is computed between the specified 'reference' array (cind) and the other arrays. Probe-specific variance estimates are robust against the choice of reference array.
#'
#' @return 
#'   fcmat: Probes x arrays preprocessed differential expression matrix.
#'   cind: Specifies which array in abatch was selected as a reference in calculating probe-level differential expression.
#'   cdf: The CDF environment used in the analysis.
#'   set.inds: Indices for probes in each probeset, corresponding to the rows of fcmat.
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords methods

RPA.preprocess <- function (abatch, 
                            bg.method = "rma",
                            normalization.method = "quantiles.robust",
                            cdf = NULL, cel.files = NULL, cel.path = NULL, quantile.basis = NULL)
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
  if (is.null(quantile.basis)) {
    abatch <- normalize(abatch, method = normalization.method)
  } else {
    # Normalize by forcing the pre-calculated quantile basis
    pm(abatch) <- set.quantiles(pm(abatch), quantile.basis)
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


determine.batch.size <- function (N, batch.size, verbose) {

  if ( is.null(batch.size) ) {
    batch.size <- min(N, 100)
  }
  if (batch.size < 3) {
    warning("Minimum batch.size is 3. Setting batch.size = 3.")
    batch.size <- 3
  }

  batch.size

}

#' get.batches
#' Split data into batches
#'
#' @param items A vector of items to be splitted into batches.
#' @param batch.size Batch size. The last batch may contain less elements than the other batches which have batch.size elements each.
#' @param shuffle Split the elements randomly in the batches.
#'
#' @return A list. Each element corresponds to one batch and contains a vector listing the elements in that batch.
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords utilities

get.batches <- function (items, batch.size = NULL, shuffle = FALSE) {

  # If batch size not given, or too small, determine automatically:	    
  batch.size <- determine.batch.size(length(items), batch.size)

  # Random ordering for the items?
  if ( shuffle ) { items <- sample(items) }

  # N elements into batches of size batch.size
  # last batch can be smaller

  if (length(items) == 1 && is.numeric(items)) {
    N <- length(items)
    items <- 1:N
  } else {
    N <- length(items)
  }

  if (N < batch.size) {
     warning("batch.size > N, setting batch.size = N.")
     batch.size <- N
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

