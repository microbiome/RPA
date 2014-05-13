library(RPA)
fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names=T); for (f in fs) {source(f)}

# Specify CEL files
cels <- list.celfiles("~/data/Kettunen13-LungCELs/Affy-data", full.names = T)
cels <- cels[-grep("HWik_z23T_110803.CEL", cels)]
cels <- cels[-grep("Hwik_z252_Barray_scan2_041203.CEL", cels)]
cels <- cels[-grep("HWik_z45N_300603.CEL", cels)]
cels <- cels[-grep("HWik_z45T_300603.CEL", cels)]
cels <- cels[-grep("HWik_z80T_310703.CEL", cels)]
cels <- cels[-grep("HWik_z80T_31scan2.CEL", cels)] 
cels <- cels[-grep("z170\\(2\\)-040203-HWik.CEL", cels)]

                cel.path = NULL;
               cel.files = cels;
                    sets = NULL;
                     cdf = NULL; 
               bg.method = "rma";                              
                  priors = list(alpha = 1, beta = 1);
                 epsilon = 1e-2;
                mc.cores = 1;
                 verbose = TRUE;                          
                 shuffle = TRUE;
              batch.size = Inf; 
                 batches = NULL; 
          quantile.basis = NULL; 
            save.batches = TRUE;
	    save.batches.dir = "."; 
	    keep.batch.files = FALSE; 
	    unique.run.identifier = NULL;
	    rseed = 231


  ###############################################################

  # FIXME: make this working also with ready-made affybatches

  # Add a unique identifier for this RPA run
  if (is.null(unique.run.identifier)) {
    unique.run.identifier <- paste("RPA-run-id-", rnorm(1), sep = "")
  }

  if (!is.null(cel.path) && is.null(cel.files)) {
    message(paste("Preprocessing all CEL files from", cel.path))
    # Randomize the order to avoid biases in batch handling
    cel.files <- sample(list.celfiles(cel.path, full.names = T))
  }

  if (is.null(batch.size)) {
    if ( verbose ) { message("Determining batch size") }
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

    # Create the output directory if necessary
    if (length(dir(save.batches.dir)) == 0) { 
      message(paste("Creating directory", save.batches.dir))
      system(paste("mkdir ", save.batches.dir)) 
    }

    message(paste("Storing intermediate batches in directory", save.batches.dir, "with the identifier", unique.run.identifier))

    blf <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-batchlist.RData", sep = "")
    if ( verbose ) { message(paste("Saving batch list into file: ", blf)) }
    save(batches, file = blf)    
  }
  
  ###############################################################

  if (is.null(quantile.basis)) {
    message("Calculating the basis for quantile normalization")    
    quantile.basis <- qnorm.basis.online(batches, bg.method, cdf, save.batches = save.batches, batch.size, verbose = verbose, save.batches.dir = save.batches.dir, unique.run.identifier = unique.run.identifier)
  }
  if (save.batches) {
    quantile.file <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-quantiles.RData", sep = "")
    if (verbose) { message(paste("Saving quantile basis into file: ", quantile.file)) }
    save(quantile.basis, file = quantile.file)
  }

  ###############################################################
source("estimate.hyper.R")    
  hyper.parameters <- estimate.hyperparameters(sets, priors,
                                                  batches, cdf,
                                                  quantile.basis,
                                                  bg.method, epsilon,
                                                  load.batches = save.batches,
                                                  save.hyperparameter.batches = save.batches,
                                                  mc.cores = 1,
                                                  verbose = TRUE, 	    
						  save.batches.dir = save.batches.dir, 
						  unique.run.identifier = unique.run.identifier)
  
  if (save.batches) {
    hyper.file <- paste(save.batches.dir, "/", unique.run.identifier, "-RPA-hyperparameters.RData", sep = "")
    if ( verbose ) { message(paste("Saving hyperparameters into file: ", hyper.file)) }
    save(hyper.parameters, file = hyper.file)
  }

  ###############################################################

  if ( verbose ) { message("Estimating affinities..") }  
  affinities <- get.affinities(sets, hyper.parameters, cdf, mc.cores, quantile.basis, bg.method, batches, load.batches = TRUE, save.batches, save.batches.dir, unique.run.identifier, verbose)

  if ( verbose ) { message("Collect hyperparameters..") } # from batch files
  hyper.parameters <- get.probe.parameters(affinities, unique.run.identifier, save.batches.dir, mode = "list") 
  hyper.parameters.evolution <- collect.hyperparameters(batches, unique.run.identifier, save.batches.dir, batch.size)

  # -----------------------------------------------------------------

  # Calculate the final summary estimates
  message("Summarizing probesets")
  emat <- summarize.batches(sets = sets, 
       	  		    variances = hyper.parameters$tau2, 
			    batches = batches, 
			    load.batches = save.batches, 
		 	    mc.cores = mc.cores, 
			    cdf = cdf, 
			    bg.method = bg.method, 
			    quantile.basis = quantile.basis, 
			    verbose = verbose, 
			    save.batches.dir = save.batches.dir, 
			    unique.run.identifier = unique.run.identifier, 
			    save.batches = save.batches, 
			    affinities = affinities)

  # -----------------------------------------------------------------

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

  ###############################################################

  # Store parameters
  params <- c(  cel.path = cel.path,
               cel.files = cel.files,
                    sets = sets,
                     cdf = cdf, 
               bg.method = bg.method,                              
                  priors = priors,
                 epsilon = epsilon,
                mc.cores = mc.cores,
                 verbose = verbose,                          
                 shuffle = shuffle,
              batch.size = batch.size, 
                 batches = batches, 
          quantile.basis = quantile.basis, 
            save.batches = save.batches,
	    save.batches.dir = save.batches.dir, 
	    keep.batch.files = keep.batch.files, 
	    unique.run.identifier = unique.run.identifier,
	    rseed = rseed)

  ###############################################################

  # Garbage collection
  gc()

eset.old <- list(expressionSet = eset, hyper.parameters = hyper.parameters, hyper.parameters.evolution = hyper.parameters.evolution, params = params, sessionInfo = sessionInfo())
