#' @title summarize.batches
#' @description Summarize batches.
#' @param sets Probesets to summarize
#' @param probe.parameters Optional probe parameters, including priors.
#' @param batches Data batches for online learning
#' @param load.batches Logical. Load precalculated data for the batches.
#' @param mc.cores Number of cores for parallel computation
#' @param cdf CDF for alternative probeset definitions
#' @param bg.method Background correction method
#' @param normalization.method Normalization method
#' @param verbose Print progress information
#' @param save.batches.dir Specify the output directory for temporary batch saves.
#' @param unique.run.identifier Define identifier for this run for naming the temporary batch files. By default, a random id is generated.
#' @param save.batches Save batches?
#' @param set.inds Probeset indices
#' @param speedup Speed up calculations with approximations.
#' @param summarize.with.affinities Use affinity estimates in probe summarization step. Default: FALSE.
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
summarize.batches <- function (sets = NULL, probe.parameters = list(), batches, load.batches = FALSE, mc.cores = 1, cdf = NULL, bg.method = "rma", normalization.method = "quantiles", verbose = TRUE, save.batches.dir = ".", unique.run.identifier = NULL, save.batches = FALSE, set.inds, speedup = FALSE, summarize.with.affinities = FALSE) {

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
    # Was in NAR version before speedup
    if (load.batches && !speedup) {
      batch.file <- paste(save.batches.dir, "/", unique.run.identifier, "-", names(batches)[[i]], ".RData", sep = "")
      if (verbose) { message(paste("Load preprocessed data for this batch from: ", batch.file)) }
      load(batch.file) # batch

      # Get probes x samples matrices for each probeset 
      # No need to remove the reference sample for d.update here 
      # in the summarization step!
      # (this returns quantile-normalized log2 data:)
      #if ( verbose ) { message("Extract probe-level data") }       
      q <- get.probe.matrix(cels = batch.cels, cdf = cdf, quantile.basis = probe.parameters$quantile.basis, bg.method = bg.method, batch = batch, verbose = verbose)

    } else if (load.batches && speedup) {
      batch.file <- paste(save.batches.dir, "/", unique.run.identifier, names(batches)[[i]], "-hyper.RData", sep = "")
      if (verbose) { message(paste("Load preprocessed data for this batch from: ", batch.file)) }
      load(batch.file) # q
    }
    
    # Summarize this batch
    res <- summarize.batch(q = q, 
    	   		   set.inds = set.inds, 	
    	   		   probe.parameters = probe.parameters,	   
			   verbose = verbose, 
			   mc.cores = mc.cores, 
			   summarize.with.affinities = summarize.with.affinities) 

    # Store the summaries
    emat[sets, batch.cels] <- res$exprs

  }

  emat

}  


