# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' probe.performance
#'
#' Provide a table of probe-level parameter estimates (affinity and stochastic noise) for RPA output.
#'
#' @param rpa.object Object of the rpa class (output from functions rpa or rpa.online)
#' @param abatch Optional: Affybatch for the rpa object, if not provided in the rpa project.
#  @param sort Sort the probes according to their stochastic noise level
#' @param sets Specify the probesets to include in the output. Default: All probesets
#' @return Data frame of probe-level parameter estimates
#' @export
#'
#' @references
#' See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # library(affydata) data(Dilution); rpa.results <- RPA.pointestimate(Dilution); tab <- probe.parameters(rpa.results); df <- df[order(abs(df$variance), decreasing = TRUE),]
#' @keywords utilities

probe.performance <- function (rpa.object, abatch = NULL, sets = NULL) {

  if (!is.null(rpa.object$abatch)) {
    # message("Using the affybatch from rpa.object")
    abatch <- rpa.object$abatch
  }

  if (is.null(sets)) {
    # Define the probesets to check
    sets <- rpa.object$sets
  }

  if (!all(sets %in% rpa.object$sets)) {
    warning("Not all sets in rpa.object, considering only the overlapping sets.")	
    sets <- intersect(sets, rpa.object$sets)
  }

  # Probe affinity effects
  af <- unlist(lapply(sets, function (set) {rpa.object[[set]]$affinity}))

  # Probe-specific noise (variance)
  s2 <- unlist(lapply(sets, function (set) {rpa.object[[set]]$tau2}))

  # PM probe indices
  pmind <- unlist(pmindex(abatch)[sets])

  # Probe effect table
  df <- data.frame(list(pmindex = pmind, affinity = af, variance = s2))

  df

}

#' probetable
#'
#' Convert probe-level hyperparameter lists into a table format.
#' 
#' @param hyper.parameters A list with alpha, betas, variances and affinities
#'
#' @export
#'
#' @references
#' See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # df <- probetable(hyper.parameters)
#' @keywords utilities

probetable <- function (hyper.parameters) {

  # Arrange probe-level parameters into a table
  sets <- unlist(sapply(1:length(hyper.parameters$mu), function (k) {rep(names(hyper.parameters$mu)[[k]], length(hyper.parameters$mu[[k]]))}))
  inds <- unlist(sapply(1:length(hyper.parameters$mu), function (k) {1:length(hyper.parameters$mu[[k]])}))
  df <- data.frame(list(probeset = sets, probe.index = inds))
  df$alpha <- rep(hyper.parameters$alpha, nrow(df)) # probe hyperparameter alpha_j (same for all probes; function of sample size T)
  df$beta <- unlist(hyper.parameters$betas)     	  # probe hyperparameter beta_j
  df$tau2 <- unlist(hyper.parameters$tau2) 	  # probe variance tau2_j (calculated from alpha, beta)
  df$mu <- unlist(hyper.parameters$mu)  	  # probe affinity mu_j   (has a prior from tau2)
  df

}


#' get.probe.parameters
#'
#' Get probe-level hyperparameter from batch files
#' 
#' @param affinities probe affinities
#' @param unique.run.identifier Batch file identifier string
#' @param save.batches.dir Batch file directory
#' @param mode "list" or "table"
#'
#' @export
#'
#' @references
#' See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # df <- get.probe.parameters(unique.run.identifier, save.batches.dir = ".", mode = "list")
#' @keywords utilities
  
get.probe.parameters <- function (affinities, unique.run.identifier, save.batches.dir = ".", mode = "list") {

  # Final estimated hyperparameters (alpha, betas, variances)
  load(paste(save.batches.dir, "/", unique.run.identifier, "-RPA-hyperparameters.RData", sep = "")) # hyper.parameters
  hyper.parameters$mu <- affinities

  if (mode == "table") {
     hyper.parameters <- probetable(hyper.parameters)
  } 

  hyper.parameters

}



#' collect.hyperparameters
#'
#' Collect probe-level parameters during online-learning from the batch files.
#' 
#' @param unique.run.identifier Batch file identifier string
#' @param save.batches.dir Batch file directory
#' @param batch.size batch size
#'
#' @export
#'
#' @references
#' See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # hpe <- collect.hyperparameters(unique.run.identifier, save.batches.dir, batch.size)
#' @keywords utilities
  
collect.hyperparameters <- function (batches, unique.run.identifier, save.batches.dir, batch.size) {

  alpha <- NULL
  betas <- NULL			

  # Hyperparameter files
  fs <- list.files(save.batches.dir, pattern = unique.run.identifier, full.names = TRUE)
  fs <- fs[grep("hyper.RData", fs)]
  o <- order(as.numeric(sapply(strsplit(unlist(strsplit(fs, "-hyper.RData")), "-"), function (x) {x[[length(x)]]})))
  fs <- fs[o]

  # Final estimated hyperparameters
  final.hp <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-hyperparameters.RData", sep = "") # hyper.parameters
  load(final.hp)

  # Hyperparameter evolution
  alphas <- c()
  betamat <- NULL
  for (f in fs) {
    load(f)
    alphas[[f]] <- alpha  
    betamat <- cbind(betamat, unlist(betas)) 
  }
  # varmat <- t(t(betamat)/alphas)

  colnames(betamat) <- names(alphas)

  # cumulative sample size from batches
  cumulative.N <- batch.size * (1:length(fs))
  cumulative.N[[length(cumulative.N)]] <- length(unlist(batches)) # alpha = N/2 + prior
  names(cumulative.N) <- names(alpha)

  # Parameter evolution
  evolution <- list()
  evolution$alpha <- alphas
  evolution$beta <- betamat
  #evolution$tau2 <- varmat
  evolution$N <- cumulative.N

  evolution

}
