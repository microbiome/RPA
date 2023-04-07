#' @title Probe performance
#' @description Provide a table of probe-level parameter estimates (affinity and stochastic noise) for RPA output.
#' @param probe.parameters List with affinities and variances for the probesets
#' @param abatch Affybatch used in the analysis
#' @param sets Specify the probesets to include in the output. Default: All probesets
#' @return Data frame of probe-level parameter estimates
#' @export
#' @references
#' See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
probe.performance <- function (probe.parameters, abatch, sets = NULL) {

  if (is.null(sets)) {
    # Define the probesets to check
    sets <- names(probe.parameters[[1]])
  }

  if (!all(sets %in% names(probe.parameters[[1]]))) {
    warning("Not all sets in rpa.object, considering only the overlapping sets.")	
    sets <- intersect(sets, names(probe.parameters[[1]]))
  }

  # Probe affinity effects
  af <- unlist(lapply(sets, function (set) {probe.parameters$affinity[[set]]}))

  # Probe-specific noise (variance)
  tau2 <- unlist(lapply(sets, function (set) {probe.parameters$tau2[[set]]}))

  # PM probe indices
  pmind <- unlist(pmindex(abatch)[sets])

  # Probe effect table
  df <- data.frame(list(pmindex = pmind, affinity = af, variance = tau2))

  df

}

#' probetable
#'
#' Convert probe-level hyperparameter lists into a table format.
#' 
#' @param probe.parameters A list with alpha, betas, variances and affinities
#'
#' @export
#'
#' @references
#' See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # df <- probetable(probe.parameters)
#' @keywords utilities

probetable <- function (probe.parameters) {

  # Arrange probe-level parameters into a table
  sets <- unlist(sapply(1:length(probe.parameters$affinity), function (k) {rep(names(probe.parameters$affinity)[[k]], length(probe.parameters$affinity[[k]]))}))
  inds <- unlist(sapply(1:length(probe.parameters$affinity), function (k) {1:length(probe.parameters$affinity[[k]])}))
  df <- data.frame(list(probeset = sets, probe.index = inds))
  df$alpha <- rep(probe.parameters$alpha, nrow(df)) # probe hyperparameter alpha_j (same for all probes; function of sample size T)
  df$betas <- unlist(probe.parameters$betas)         # probe hyperparameter beta_j
  df$tau2 <- unlist(probe.parameters$tau2) 	    # probe variance tau2_j (calculated from alpha, beta)
  df$affinity <- unlist(probe.parameters$affinity)  # probe affinity mu_j   (has a prior from tau2)
  df$quantile.basis <- probe.parameters$quantile.basis
  df

}

#' probe.parameters.tolist
#'
#' Convert probe parameter table into a list format 
#' 
#' @param probe.parameters A data.frame with alpha, betas, tau2, affinities, quantile.basis
#'
#' @export
#'
#' @references
#' See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # df <- probe.parameters.tolist(probe.parameters.table)
#' @keywords utilities

probe.parameters.tolist <- function (probe.parameters) {

  if (is.data.frame(probe.parameters)) {

    message("Converting probe.parameters into list format")
    splitted <- split(probe.parameters, probe.parameters$probeset)

    probe.parameters.list <- list()    
    if (length(unique(probe.parameters$alpha)) > 1) {stop("Assuming alpha is identical everywhere.")}
    probe.parameters.list$alpha <- unique(probe.parameters$alpha)
    probe.parameters.list$betas <- lapply(splitted, function (set.data) {set.data[order(set.data$probe.index), "beta"]})
    probe.parameters.list$tau2 <- lapply(splitted, function (set.data) {set.data[order(set.data$probe.index), "tau2"]})
    probe.parameters.list$affinity <- lapply(splitted, function (set.data) {set.data[order(set.data$probe.index), "affinity"]})
    probe.parameters.list$quantile.basis <- unlist(lapply(splitted, function (set.data) {set.data[order(set.data$probe.index), "quantile.basis"]}))

  } else {
    probe.parameters.list <- probe.parameters
  }

  probe.parameters.list

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
  load(paste(save.batches.dir, "/", unique.run.identifier, "-RPA-hyperparameters.RData", sep = "")) # probe.parameters
  probe.parameters$mu <- affinities

  if (mode == "table") {
     probe.parameters <- probetable(probe.parameters)
  } 

  probe.parameters

}



#' collect.hyperparameters
#'
#' Collect probe-level parameters during online-learning from the batch files.
#' 
#' @param batches batch list
#' @param unique.run.identifier Batch file identifier string
#' @param save.batches.dir Batch file directory
#' @param save.batches Logical. Determines whether batches are available.
#' @param verbose verbose
#'
#' @export
#'
#' @references
#' See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # hpe <- collect.hyperparameters(batches, unique.run.identifier, save.batches.dir, save.batches)
#' @keywords utilities
  
collect.hyperparameters <- function (batches, unique.run.identifier, save.batches.dir, save.batches, verbose = TRUE) {

  if ( verbose ) { message("Collecting hyperparameters..") } # from batch files

  evolution <- NULL

  if (is.null(save.batches)) {

    alpha <- NULL
    betas <- NULL			

    # Hyperparameter files
    fs <- list.files(save.batches.dir, pattern = unique.run.identifier, full.names = TRUE)
    fs <- fs[grep("hyper.RData", fs)]
    o <- order(as.numeric(sapply(strsplit(unlist(strsplit(fs, "-hyper.RData")), "-"), function (x) {x[[length(x)]]})))
    fs <- fs[o]

    # Final estimated hyperparameters
  final.hp <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-hyperparameters.RData", sep = "") # probe.parameters
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
    batch.size <- length(batches[[1]])
    cumulative.N <- batch.size * (1:length(fs))
    cumulative.N[[length(cumulative.N)]] <- length(unlist(batches)) # alpha = N/2 + prior
    names(cumulative.N) <- names(alpha)

    # Parameter evolution
    evolution <- list()
    evolution$alpha <- alphas
    evolution$beta <- betamat
    #evolution$tau2 <- varmat
    evolution$N <- cumulative.N

  } 

  evolution
}
