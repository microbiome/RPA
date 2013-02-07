# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

rpa.impute <- function (dat) {
  if (sum(is.na(dat)) > 0) {
    warning(paste("Data has ", mean(is.na(dat)), " fraction of missing values: imputing"))
    dat <- t(apply(dat, 1, function (x) { y <- x; y[is.na(y)] <- rnorm(sum(is.na(y)), mean(y, na.rm = T), sd(y, na.rm = T)); y}))
  }

  dat
}


set.alpha <- function (alpha = NULL, tau2.method, P){ 
  
  # set uninformative prior if not given
  if ((tau2.method == "mean" || tau2.method == "online" || tau2.method == "robust")) {
    if (is.null(alpha)) { 
      alpha <- 1 + 1e-6
    } else if (any(alpha <= 1)) {
      stop(paste("Set alpha > 1!"))
    }
  } else if (is.null(alpha)) { 
    alpha <- 1e-6 
  }

  alpha
}

set.beta <- function (beta, tau2.method, P) {

  # if beta is scalar, set identical prior for all probes with this value
  if (is.null(beta) && !tau2.method == "robust") {
    beta <- rep.int(1e-6, P)
  } else if (is.null(beta) && tau2.method == "robust") {
    beta <- rep.int(1, P)
  } else if (length(beta) == 1) {
    beta <- rep.int(beta, P)
  } else {}

  beta

}

centerData <- function (X,rm.na = FALSE, meanvalue = NULL) {


  if (!rm.na) {
    xcenter <- colMeans(X)
    X2 <- X - rep(xcenter, rep.int(nrow(X), ncol(X)))
  } else {	
    X2 <- array(NA, dim = c(nrow(X), ncol(X)), dimnames = dimnames(X))
    for (i in 1:ncol(X)) {
      x <- X[,i]
      nainds <- is.na(x)
      xmean <- mean(x[!nainds])
      X2[!nainds,i] <- x[!nainds] - xmean 	
    }
    dimnames(X2) <- dimnames(X)
  }

  if (!is.null(meanvalue)) {
    # Shift the data so that mean gets a specified value
    X2 <- X2 + meanvalue
  }


  
  X2
}

