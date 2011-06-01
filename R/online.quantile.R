# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


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


  cat("Calculating quantile basis..\n")
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
  # replace smallest value with the smallest in quantile.basis etc.
  apply(mat, 2, function(x) { 
    x[order(x)] <- sort(quantile.basis)
    x
  })
}
