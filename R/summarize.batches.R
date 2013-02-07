# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' summarize.batches
#'
#' @param sets Probesets to summarize
#' @param variances Precalculated probe-specific variances
#' @param batches Data batches for online learning
#' @param load.batches Logical. Load precalculated data for the batches.
#' @param mc.cores Number of cores for parallel computation
#' @param cdf CDF for alternative probeset definitions
#' @param bg.method Background correction method
#' @param normalization.method Normalization method
#' @param verbose Print progress information
#' @param quantile.basis Basis for quantile normalization
#' @param save.batches.dir Specify the output directory for temporary batch saves.
#' @param unique.run.identifier Define identifier for this run for naming the temporary batch files. By default, a random id is generated.
#' @param save.batches Save batches?
#' @param affinities probe-specific affinities
#' @param set.inds Probeset indices
#'
#' @details Sweeps through the batches. Summarizes the probesets within each batch based on the precalculated model parameter point estimates.
#'
#' @return Expression matrix: probesets x samples.
#'
#' @export
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords utilities

summarize.batches <- function (sets = NULL, variances, batches, load.batches = FALSE, mc.cores = 1, cdf = NULL, bg.method = "rma", normalization.method = "quantiles", verbose = TRUE, quantile.basis, save.batches.dir = ".", unique.run.identifier = NULL, save.batches = FALSE, affinities, set.inds) {

  if ( verbose ) { message("Pick PM indices") }
  if (is.null(sets)) {sets <- names(set.inds)}
  
  # Initialize expression matrix
  cel.files <- unlist(batches)
  emat <- array(NA, dim = c(length(sets), length(cel.files)))
  rownames(emat) <- sets
  colnames(emat) <- cel.files # sapply(strsplit(cel.files, "/"), function (x) {x[[length(x)]]})

  # Hyperparameters for probe variances have been estimated
  # Get probeset-level signal estimate by online sweep with the fixed variances  
  for (i in 1:length(batches)) {

    if (verbose) {message(paste("Summarizing batch", i, "/", length(batches)))}

    # CEL files for this batch
    batch.cels <- batches[[i]] 

    # Get background corrected, quantile normalized, and logged probe-level matrix
    batch <- NULL
    if (load.batches) {
      batch.file <- paste(save.batches.dir, "/", unique.run.identifier, "-", names(batches)[[i]], ".RData", sep = "")
      if (verbose) { message(paste("Load preprocessed data for this batch from: ", batch.file)) }
      load(batch.file) # batch
    }
    
    # Get probes x samples matrices for each probeset 
    # No need to remove the reference sample for d.update here 
    # in the summarization step!
    # (this returns quantile-normalized log2 data:)
    if ( verbose ) { message("Extract probe-level data") }       
    q <- get.probe.matrix(cels = batch.cels, cdf = cdf, quantile.basis = quantile.basis, bg.method = bg.method, batch = batch, verbose = verbose)

    probe.parameters <- list(affinity = affinities, variance = variances)

    # Summarize this batch
    mat <- summarize.batch(q = q, set.inds = set.inds, 	
    	   		   probe.parameters = probe.parameters,	   
			   verbose = verbose, 
			   mc.cores = mc.cores) 

    # Store the summaries
    emat[sets, batch.cels] <- mat

  }

  emat

}  


