# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2011-2013 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


qnorm.basis.online <- function (cel.files, bg.method = "rma", cdf = NULL, save.batches = FALSE, batch.size = 3, verbose = TRUE, save.batches.dir = ".", unique.run.identifier = NULL) {

  message("Calculating the basis for quantile normalization")    

  # Estimate basis for quantile normalization 

  # FIXME: to speed up, add option to calculate only based on subset of data (as
  # in online.quantile function; merge these two)
  
  # Split CEL file list into batches
  if (length(cel.files[[1]]) == 1) {

    # Assume cel.files is a character vector listing CEL files
    # Create batches
    batches <- get.batches(cel.files, batch.size, shuffle = TRUE)

  } else if (length(cel.files[[1]]) > 1) {

    # Assume that cel.files is already a list of batches
    batches <- cel.files
    cel.files <- unlist(batches)
  }

  if (is.null(names(batches))) {
    names(batches) <- paste("batch", 1:length(batches), sep = "-")
  }
  
  if (save.batches) {
    message(paste("Saving background corrected data for quantile normalization into temporary files to speed up preprocessing in later steps.", sep = ""))
  }

  # Process each batch separately
  qs <- NULL

  for (i in 1:length(batches)) {

    message(paste("Calculating quantiles on batch", i, "/", length(batches)))

    abatch <- ReadAffy(filenames = batches[[i]], compress=getOption("BioC")$affy$compress.cel) 

    if (!is.null(cdf)) {
      if (verbose) {message("Setting alternative CDF environmen")}
      abatch@cdfName <- cdf
    }    

    if (verbose) {message("Background correcting...")}
    abatch <- bg.correct(abatch, bg.method, destructive = TRUE)    
    if (verbose) {message("Done.")}    
    pma <- pm(abatch)

    # Store the necessary PM probe information used by quantile normalization
    # This can speed up preprocessing considerably
    if (save.batches) {
      if (verbose) { message("Saving the batch") }
      batch <- apply(pma, 2, rank)
      colnames(batch) <- batches[[i]]
      bf <- paste(save.batches.dir, "/", unique.run.identifier, "-", names(batches)[[i]], ".RData", sep = "")
      save(batch, file = bf)
      if ( verbose ) { message(paste("Stored batch data in: ", bf)) }      
    }

    # Add to overall probe-sum for quantile estimation
    if ( verbose ) { message("Updating the quantiles...") }    
    if (is.null(qs)) {
      # Create new quantile basis
      qs <- rowSums(apply(pma, 2, sort))
    } else {
      # Add to existing quantile basis
      qs <- qs + rowSums(apply(pma, 2, sort))
    }
    if (verbose) {message("Done.")}    
  }

  # Quantile basis is average (sum/n) over the individual arrays
  # Finally, the data is presented in log2
  #if (verbose) {message("")}
  basis <- log2(qs/length(cel.files))

}



#' online.quantile
#' Quantile normalization tools for online preprocessing. Estimate quantiles for quantile normalization based on subset of the data (random, or specified by the user).
#'
#' @param abatch AffyBatch
#' @param n Numeric: number of random samples to use to define quantile basis. Vector: specify samples to be used in quantile basis calculation.
#'
#' @details "online.quantile": Ordinary quantile normalization is exhaustively memory-consuming in alrge data sets. Then the quantiles can be calculated based on subset of the data to allow efficient normalization. This function can also be used to investigate effect of subset size to convergence of the quantile estimates;"qnorm.basis.online": sweeps through the data in batches to calculate the basis for quantile normalization (average over sorted profiles).
#'
#' @return "online.quantile": AffyBatch; "qnorm.basis.online": a vector containing the basis for quantile normalization.
#'
#' @export
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # 
#' @keywords utilities

online.quantile <- function (abatch, n) {

  # Pick random subset of abatch and learn quantile normalization
  # basis from this. Then apply the same basis for the other (new)
  # arrays. Test how many samples are sufficient for basis
  # convergence.

  if (length(n) == 1) {
    if (n == 1) {stop("Provide n > 1 for quantile basis estimation!")}
    # Pick random sample of n arrays
    s <- sample(ncol(exprs(abatch)), n)
  } else {
    # If n is a vector, select the samples indicated by n
    s <- n
  }

  message("Calculating quantile basis..\n")
  # FIXME: provide robust version
  qb <- get.quantile.basis(pm(abatch)[,s])

  cat("Replace data with quantiled values..\n")
  pm(abatch) <- set.quantiles(pm(abatch), qb)

  abatch

}  

get.quantile.basis <- function (mat) {
  sort(rowMeans(apply(mat, 2, sort)))
}


set.quantiles <- function (mat, quantile.basis) {
  apply(mat, 2, function(x) {quantile.basis[rank(x)]})
}




