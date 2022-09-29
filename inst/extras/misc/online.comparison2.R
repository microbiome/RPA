library(RPA2)
fs <- list.files("RPA/pkg/R/", full.names=T); for (f in fs) {source(f)}

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
	    save.batches.dir = "."; 
	    keep.batch.files = FALSE; 
	    unique.run.identifier = paste("RPA-run-id-", rnorm(1), sep = "");
	    rseed = 23

  # Always store the batches temporarily
  save.batches <- TRUE

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

  # Optionally save batch list
  blf <- saving.batchlist(batches, save.batches, save.batches.dir, unique.run.identifier, verbose)
  
  # ---------------------------------------------------------------------------

  # CALCULATE THE BASIS FOR QUANTILE NORMALIZATION
  if (is.null(quantile.basis)) {
    quantile.basis <- qnorm.basis.online(batches, bg.method, cdf, save.batches = save.batches, verbose = verbose, save.batches.dir = save.batches.dir, unique.run.identifier = unique.run.identifier)
  }

  # Optionally save quantile basis
  qf <- saving.quantile.basis(quantile.basis, save.batches, save.batches.dir, unique.run.identifier, verbose)

  # ---------------------------------------------------------------------------------

  # HYPERPARAMETER ESTIMATION

  if (verbose) { message("Get probeset indices") }
  set.inds <- get.set.inds(batches[[1]][1:2], cdf, sets)

  hyper.parameters <- estimate.hyperparameters(sets, priors,
                                                  batches, cdf,
                                                  quantile.basis,
                                                  bg.method, epsilon,
                                                  load.batches = save.batches,
                                                  save.hyperparameter.batches = save.batches,
                                                  mc.cores = mc.cores,
                                                  verbose = verbose, 	    
						  save.batches.dir = save.batches.dir, 
						  unique.run.identifier = unique.run.identifier, 
						  set.inds = set.inds)

stop("HERE")

  hf <- saving.hyperparameters(hyper.parameters, save.batches, save.batches.dir, unique.run.identifier, verbose)

  # -----------------------------------------------------------------------------------

  # AFFINITY ESTIMATION
  affinities <- get.affinities(sets, hyper.parameters, cdf, mc.cores, quantile.basis, bg.method, batches, load.batches = save.batches, save.batches, save.batches.dir, unique.run.identifier, verbose, set.inds)

  # Add affinities into hyperparameter list
  hyper.parameters$mu <- affinities

  # ----------------------------------------------------------------------------------------------
  
  # COLLECTING HYPERPARAMETERS ACROSS ALL BATCHES TO INVESTIGATE HYPERPARAMETER EVOLUTION
  hyper.parameters.evolution <- collect.hyperparameters(batches, unique.run.identifier, save.batches.dir, save.batches)

  # -----------------------------------------------------------------------------------

  # PROBESET SUMMARIZATION

  # Given the affinities, calculate the final summary estimates
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
			    affinities = affinities, set.inds = set.inds)

  # Arrange CEL files in the original order and Coerce expression
  # values in the rpa object into an ExpressionSet object and return
  # expression set object
  eset <- new("ExpressionSet", assayData = list(exprs = emat[, cel.files])) 

  # ---------------------------------------------------------------------------

  # Remove temporary batch files and run garbage collection
  tmp <- batch.cleanup(keep.batch.files, unique.run.identifier, save.batches.dir)

  # -------------------------------------------------------------------------------------

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
	      batch.size = length(batches[[1]]),
                 batches = batches, 
          quantile.basis = quantile.basis, 
            save.batches = save.batches,
        save.batches.dir = save.batches.dir, 
	keep.batch.files = keep.batch.files, 
   unique.run.identifier = unique.run.identifier,
	           rseed = rseed, 
	     sessionInfo = sessionInfo())

eset.new <- list(expressionSet = eset, hyper.parameters = hyper.parameters, hyper.parameters.evolution = hyper.parameters.evolution, params = params)
