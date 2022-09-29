library(RPA)
fs <- list.files("~/Rpackages/RPA/github/RPA/pkg/R/", full.names = T)
for (f in fs) {source(f)}
cels <- list.celfiles("CEL", full.names = T)


                cel.path = NULL;
               cel.files = cels; #sample(cels, 400);
                    sets = NULL;
                     cdf = NULL; 
               bg.method = "rma";                              
        probe.parameters = list(alpha = 1, beta = 1);
                 epsilon = 1e-2;
                mc.cores = 4;
                 verbose = TRUE;                          
                 shuffle = TRUE;
              batch.size = 200; 
                 batches = NULL; 
	save.batches.dir = "."; 
	keep.batch.files = FALSE; 
   unique.run.identifier = "TESTING";
	           rseed = 23; 
		 speedup = FALSE; 
		 summarize.with.affinities = FALSE
		 
  # Get CEL file list
  if (!is.null(cel.path) && is.null(cel.files)) {
    message(paste("Preprocessing all CEL files from", cel.path))
    cel.files <- list.celfiles(cel.path, full.names = T)
  }

  # Define batches
  if ( is.null(batches) ) {
    message("Split CEL file list into batches")
    batches <- get.batches(cel.files, batch.size, shuffle)
  }

  # Optionally save batch 
  if (!speedup) {
    blf <- saving.batchlist(batches, save.batches.dir, unique.run.identifier, verbose)
  }

  # ---------------------------------------------------------------------

  # CALCULATE THE BASIS FOR QUANTILE NORMALIZATION

  # TODO can be sped up by calculating quantiles based on data subset only
  # and then calculating bg corrections as part of hyperparameter estimation
  # then we will have one less round of save/load iterations which 
  # will save considerable time with large data sets
  if (is.null(probe.parameters$quantile.basis)) {
    if (speedup) {
      nq <- 500
      message(paste("Estimating quantiles based on random subset of", nq, "CEL files to speed up preprocessing.."))
      nbs <- min(nq, length(cel.files))
      qbatches <- get.batches(sample(cel.files, nbs), batch.size, shuffle)
      #qbatches <- list(qbatch = sample(cel.files, nbs))
      probe.parameters$quantile.basis <- qnorm.basis.online(qbatches, bg.method, cdf, save.batches = FALSE, verbose = verbose)
      
    } else {
      probe.parameters$quantile.basis <- qnorm.basis.online(batches, bg.method, cdf, save.batches = TRUE, verbose = verbose, save.batches.dir = save.batches.dir, unique.run.identifier = unique.run.identifier)
      
    }
  }

  # Optionally save quantile basis
  if (!speedup) {    
    quantile.basis <- probe.parameters$quantile.basis
    qf <- saving.quantile.basis(quantile.basis, save.batches.dir, unique.run.identifier, verbose)
  }

  # ---------------------------------------------------------------------

  # HYPERPARAMETER ESTIMATION
  if (verbose) { message("Get probeset indices") }
  set.inds <- get.set.inds(batches[[1]][1:2], cdf, sets)

  message("Estimating probe.parameters")

    if (speedup) {
      probe.parameters <- estimate.hyperparameters(sets = sets, 
      		       	  		     probe.parameters = probe.parameters,
                                              batches = batches, 
					      cdf = cdf,
                                              bg.method = bg.method, 
					      epsilon = epsilon,
                                              load.batches = FALSE,
                                        save.hyperparameter.batches = TRUE,
                                                  mc.cores = mc.cores,
                                                  verbose = verbose, 	    
					save.batches.dir = save.batches.dir, 
				unique.run.identifier = unique.run.identifier, 
						  set.inds = set.inds)

    } else {
      probe.parameters <- estimate.hyperparameters(sets, probe.parameters,
                                              batches, cdf,
                                              bg.method, epsilon,
                                              load.batches = TRUE,
                                        save.hyperparameter.batches = TRUE,
                                                  mc.cores = mc.cores,
                                                  verbose = verbose, 	    
					save.batches.dir = save.batches.dir, 
				unique.run.identifier = unique.run.identifier, 
						  set.inds = set.inds)

    }


  hf <- saving.hyperparameters(probe.parameters, save.batches.dir, 
     			       unique.run.identifier, verbose)

  # ----------------------------------------------------------------

  # AFFINITY ESTIMATION
  if (is.null(probe.parameters$affinity)) {
    #if (speedup) {
    #  message("Skipping affinity estimation")
    #  probe.parameters$affinity <- NULL
    #} else {
      message("Estimating affinities")
      probe.parameters$affinity <- get.affinities(sets, probe.parameters, cdf, mc.cores, bg.method, batches, load.batches = TRUE, save.batches = TRUE, save.batches.dir = save.batches.dir, unique.run.identifier = unique.run.identifier, verbose = verbose, set.inds = set.inds)
    #}
  } else {
    warning("Affinities provided in rpa.online function call through probe.parameters. Using these precalculated affinities and skipping affinity estimation!") 
  }

  # ------------------------------------------------------------------
  
  # COLLECTING HYPERPARAMETERS ACROSS ALL BATCHES TO INVESTIGATE 
  # HYPERPARAMETER EVOLUTION
  probe.parameters.evolution <- NULL
  if (speedup) {
    message("Skipping hyperparameter converence analysis to speed up...")
  } else {
    message("Investigate parameter convergence...")
    probe.parameters.evolution <- collect.hyperparameters(batches, unique.run.identifier, save.batches.dir, save.batches = TRUE)
    message("Parameter convergence list done.")
  }

  # --------------------------------------------------------------

  # PROBESET SUMMARIZATION

  # Given the affinities, calculate the final summary estimates
  message("Summarizing probesets")
  emat <- summarize.batches(sets = sets, 
       	  		    probe.parameters = probe.parameters,
			    batches = batches, 
			    load.batches = TRUE, 
		 	    mc.cores = mc.cores, 
			    cdf = cdf, 
			    bg.method = bg.method, 
			    verbose = verbose, 
			    save.batches.dir = save.batches.dir, 
			    unique.run.identifier = unique.run.identifier, 
			    save.batches = TRUE, 
			    set.inds = set.inds, 
			    speedup = speedup, 
			    summarize.with.affinities = summarize.with.affinities)

  # Arrange CEL files in the original order and Coerce expression
  # values in the rpa object into an ExpressionSet object and return
  # expression set object
  eset <- new("ExpressionSet", assayData = list(exprs = emat[, cel.files])) 

  # ------------------------------------------------------------------

  # Remove temporary batch files and run garbage collection
  tmp <- batch.cleanup(keep.batch.files, unique.run.identifier, save.batches.dir)

  # -------------------------------------------------------------------

  message("Save parameters")
  params <- c(  cel.path = cel.path,
               cel.files = cel.files,
                    sets = sets,
                     cdf = cdf, 
               bg.method = bg.method,                              
                 epsilon = epsilon,
                mc.cores = mc.cores,
                 verbose = verbose,                          
                 shuffle = shuffle,
	      batch.size = length(batches[[1]]),
                 batches = batches, 
        save.batches.dir = save.batches.dir, 
	keep.batch.files = keep.batch.files, 
   unique.run.identifier = unique.run.identifier,
	           rseed = rseed, 
		   speedup = speedup,
 summarize.with.affinities = summarize.with.affinities,
	     sessionInfo = sessionInfo())


res <-  list(expressionSet = eset, probe.parameters = probe.parameters, probe.parameters.evolution = probe.parameters.evolution, params = params)


