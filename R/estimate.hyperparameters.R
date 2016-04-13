#' @title estimate.hyperparameters
#' @description Hyperparameter estimation.
#' @param sets Probesets to handle. All probesets by default.
#' @param probe.parameters User-defined priors. May also include quantile.basis
#' @param batches Data batches for online learning
#' @param cdf CDF probeset definition file
#' @param bg.method Background correction method
#' @param epsilon Convergence parameter
#' @param load.batches Logical. Load preprocessed data whose identifiers are picked from names(batches). Assuming that the same batch list (batches) was used to create the files in online.quantiles function.
#' @param save.hyperparameter.batches Save hyperparameters for each batch into files using the identifiers with batch name with -hyper.RData suffix.
#' @param mc.cores Number of cores for parallel computation
#' @param verbose Print progress information
#' @param normalization.method Normalization method
#' @param save.batches.dir Specify the output directory for temporary batch saves.
#' @param unique.run.identifier Define identifier for this run for naming the temporary batch files. By default, a random id is generated.
#' @param set.inds Probeset indices
#'
#'@return alpha: Hyperparameter alpha (same for all probesets); betas: Hyperparameter beta (probe-specific); variances: Probe-specific variances (beta/alpha)
#'
#' @references See citation("RPA") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples # 
#' @keywords utilities

estimate.hyperparameters <- function (sets = NULL, 
			  probe.parameters = list(alpha = 2, beta = 1), 
			    	      batches, cdf = NULL, 
				      bg.method = "rma", 
				      epsilon = 1e-2,
                                      load.batches = FALSE,
                                      save.hyperparameter.batches = FALSE, 
				      mc.cores = 1, 
				      verbose = TRUE, 
				      normalization.method = "quantiles",
				      save.batches.dir = ".", 
				      unique.run.identifier = NULL,
				      set.inds = set.inds)
{

  # Hyperparameter estimation through batches  

  if ( is.null(sets) ) { sets <- names(set.inds) }
  
  # Initialize hyperparameters
  # Note: alpha is scalar and same for all probesets 
  # alpha <- alpha + N/2 at each batch
  if (verbose) { message("Initialize priors") }
  alpha <- probe.parameters$alpha # initialize 
  if (length(unlist(probe.parameters$beta)) == 1) { # If input beta is a scalar
    betas <- mclapply(sets, function (set) { 
      rep(probe.parameters$beta, length(set.inds[[set]])) 
    }, mc.cores = mc.cores)
    names(betas) <- sets
    probe.parameters$beta <- betas  
  }
  
  for (i in 1:length(batches)) {

    if (verbose) {
      message(paste("Updating hyperparameters; batch", i, "/", length(batches)))
    }
    
    # Load batch with bgc, ordered data
    batch <- NULL
    if (load.batches) {
      batch.file <- paste(save.batches.dir, "/", unique.run.identifier, "-", names(batches)[[i]], ".RData", sep = "")      
      if (verbose) { message(paste("Load batch from file:", batch.file)) }
      load(batch.file) # batch
    }

    # Pick quantile.basis from probe parameters if given
    quantile.basis <- probe.parameters$quantile.basis

    # Get background corrected, quantile normalized, log2 probe-level matrix
    if ( verbose ) { message("Pick probe-level values") }    
    q <- get.probe.matrix(cels = batches[[i]], cdf = cdf, 
      	 		       quantile.basis = quantile.basis,
                               bg.method = bg.method, 
			       normalization.method = normalization.method, 
			       batch = batch,
                               verbose = verbose)

    # -----------------------------------------------
    
    # Update hyperparameters
    if (is.null(probe.parameters$tau2)) {
      hp <- updating.hyperparameters(q, set.inds, verbose, mc.cores = mc.cores, alpha, betas, epsilon)
      alpha <- hp$alpha
      betas <- hp$betas
      tau2  <- hp$s2s
    } else {

      warning("Probe variances already provided in the input argument probe.parameters. Skipping hyperparameter estimation and using the predefined variances!")

      alpha <- probe.parameters$alpha
      betas <- probe.parameters$betas
      tau2  <- probe.parameters$tau2
    }

    bf <- saving.hyperparameter.batches(alpha, betas, save.hyperparameter.batches, save.batches.dir, unique.run.identifier, names(batches)[[i]], verbose, q)
        
  }

  # Get final estimated variances for each probeset based on hyperparameter posteriors
  # variances <- mclapply(betas, function (beta) {beta/alpha}, mc.cores = mc.cores)

  list(alpha = alpha, betas = betas, tau2 = tau2, quantile.basis = quantile.basis)  

}


saving.hyperparameter.batches <- function (alpha, betas, save.hyperparameter.batches, save.batches.dir, unique.run.identifier, nam, verbose, q) {

    if (save.hyperparameter.batches) {
      batch.file <- paste(save.batches.dir, "/", unique.run.identifier, nam, "-hyper.RData", sep = "")
      if ( verbose ) { message(paste("Save hyperparameters into file:", batch.file)) }
      save(q, alpha, betas, file = batch.file)
    }
}


