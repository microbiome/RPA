#rpa.online <- function (
#)
#{

#cel.path = NULL; cel.files = cels; sets = sets; cdf = cdf; bg.method = "rma"; priors = list(alpha = 1, beta = 1); epsilon = 1e-2; mc.cores = 4; verbose = TRUE; shuffle = TRUE; batch.size = 18; batches = NULL; quantile.basis = NULL; save.batches = TRUE; save.batches.dir = "."

  ###############################################################

  warning("rpa.online is an experimental version")

  # Add a unique identifier for this RPA run
  if (is.null(unique.run.identifier)) {
    unique.run.identifier <- paste("RPA-run-id-", rnorm(1), "-", sep = "")
  }

  
  if (save.batches) {

    # Create the output directory if necessary
    if (length(dir(save.batches.dir)) == 0) { 
      system(paste("mkdir ", save.es.dir)) 
    }

    message(paste("Storing intermediate batches in directory", save.batches.dir, "with the identifier", unique.run.identifier))

  }

  if (is.null(batch.size)) {
    if (verbose) {message("Determining batch size")}
    batch.size <- min(length(unlist(cel.files)), 100)
  }
  
  if (batch.size < 3) {
    warning("Minimum batch.size is 3. Setting batch.size = 3.")
    batch.size <- 3
  }

  if (is.null(batches)) {
    message("Split CEL file list into batches")
    batches <- get.batches(cel.files, batch.size, shuffle)
  }
  if (save.batches) {
    blf <- paste(save.batches.dir, "/", unique.run.identifier, "RPA-batchlist.RData", sep = "")
    if (verbose) {message(paste("Saving batch list into file: ", blf))}
    save(batches, file = blf)    
  }
  
  ###############################################################

  if (is.null(quantile.basis)) {
    message("Calculating the basis for quantile normalization")    
    quantile.basis <- qnorm.basis.online(batches, bg.method, cdf, save.batches = save.batches, batch.size, verbose = verbose, save.batches.dir = save.batches.dir, unique.run.identifier = unique.run.identifier)
  }
  if (save.batches) {
    quantile.file <- paste(save.batches.dir, "/", unique.run.identifier, "RPA-quantiles.RData", sep = "")
    if (verbose) { message(paste("Saving quantile basis into file: ", quantile.file)) }
    save(quantile.basis, file = quantile.file)
  }

  ###############################################################
    
  hyper.parameters <- estimate.hyperparameters(sets, priors,
                                                  batches, cdf,
                                                  quantile.basis,
                                                  bg.method, epsilon,
                                                  load.batches = save.batches,
                                                  save.hyperparameter.batches = save.batches,                                                  
                                                  mc.cores = mc.cores,
                                                  verbose = verbose, 	    
						  save.batches.dir = save.batches.dir, 
						  unique.run.identifier = unique.run.identifier)
  
  if (save.batches) {
    hyper.file <- paste(save.batches.dir, "/", unique.run.identifier, "RPA-hyperparameters.RData", sep = "")
    if (verbose) {message(paste("Saving hyperparameters into file: ", hyper.file))}
    save(hyper.parameters, file = hyper.file)
  }

  ###############################################################

  # Final ExpressioSet object 
  message("Summarizing probesets")

  emat <- summarize.batches(sets = sets, 
       	  		    variances = hyper.parameters$variances, 
			    batches = batches, 
			    load.batches = save.batches, 
			    mc.cores = mc.cores, 
			    cdf = cdf, 
			    bg.method = bg.method, 
			    quantile.basis = quantile.basis, 
			    verbose = verbose, 
			    save.batches.dir = save.batches.dir, 
			    unique.run.identifier = unique.run.identifier)
  
  ##################################################################

  if (!keep.batch.files) {
    message(paste("Removing the temporary batch files from directory", save.batches.dir, "with the identifier", unique.run.identifier))
    system(paste("rm ", save.batches.dir, "/", unique.run.identifier, "*", sep = ""))
  } else {
    message(paste("Keeping the temporary batch files in directory", save.batches.dir, "with the identifier", unique.run.identifier))
  }
  
  # Arrange CEL files in the original order and Coerce expression
  # values in the rpa object into an ExpressionSet object and return
  # expression set
  eset <- new("ExpressionSet", assayData = list(exprs = emat[, cel.files])) 

