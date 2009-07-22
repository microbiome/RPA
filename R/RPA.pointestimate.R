RPA.pointestimate <-
function (abatch, sets=NULL,
		myseed = 101, 
		priors = NULL,
		epsilon = 10^(-2), 
		cind = FALSE, sigma2.method = "var", d.method = "fast", verbose = TRUE) 
{

#
# Find posterior mode for RPA model parameters d (mean) and sigma2 (variances)
#
# INPUT: 
#
#	abatch : affybatch object, as obtained with ReadAffy()
#
#	sets : Character vector specifying probesets for which 
#	       RPA will be computed. Default: all probe sets.
#
#	myseed : random seed to be used
#
# 	priors : priors for the model parameters.
#		 an instance of the class 'rpa.priors'
#		 alpha and beta are priors for inverse gamma distribution i.e. for sigma2
#		      (uninformative as alpha, beta -> 0)
#		      NOTE: default 1e-6 used if prior has not been specified
#		      priors used only when sigma2.method = "basic"
#		 prior for d has not been implemented
#
#	epsilon: convergence threshold 
#		(maximal accepted change in all parameters)
#
#	cind : index of the control array (in the CEL file list). 
#	       If FALSE (default), control is chosen randomly.
#	       Default is FALSE; the method has appeared to be very 
#	       robust for the arbitrary choice of cind.
#
#	sigma2.method: "var", "basic"  
#
#	d.method: "fast", "basic"
#
#	verbose: print progress information. Default: TRUE.
#
#
#
# OUTPUT:
#
#	res : output list containing the following entries
#		res$d : differential gene expression estimate for each probe set 
#			between the randomly selected 'control' array and the remaining arrays
#		res$sigma2 : probe-specific variances for the probes in each probe set
#		res$cind : index of the control array in the input data
#
#
#
# EXAMPLE:
#
#	rpa.results <- RPA.pointestimate(abatch)
#


	#################################################################

	# PREPROCESSING

	#################################################################

	#Set random seed
	set.seed(myseed)

	# Preprocessing
	preproc <- RPA.preprocess(abatch, cind)

	# Pick the necessary objects
	fcmat <- preproc$fcmat # probe-level differential expressions
	cind <- preproc$cind   # index of the control array

	#################################################################

	# ESTIMATE PROBE RELIABILITY AND DIFFERENTIAL GENE EXPRESSION

	#################################################################

	# Number of arrays except control
	T <- ncol(fcmat)

	# Check names and number for the investigated probesets
	# if my.sets not specified, take all sets in abatch
	if (length(sets)==0) {sets<-geneNames(abatch)} 
	Nsets<-length(sets) 

	## Matrices for storing the results
	d.results<-matrix(NA,nr=Nsets,nc=T)
	sigma2.results<-vector(length=Nsets,mode="list")

	names(sigma2.results)<-sets
	rownames(d.results)<-sets
	colnames(d.results)<-colnames(fcmat)

	for (i in seq(Nsets)) {
	
		if (verbose) {
			print(paste("Computing probeset",sets[[i]],":",i,"/",Nsets,"..."))
		}

		#find probe (pm) indices for this set
		pmindices<-pmindex(abatch,sets[[i]])[[1]]
	
		#Number of probes for this probeset
		P<-length(pmindices)

		# Set initial guesses for sigma (P-vector of probe-wise stds)
		sigma2.guess<-rep(1,P) #give equals for all probes initially
		
		#Get chips x probes matrix of probe-wise fold-changes
		S<-t(fcmat[pmindices,])
		
		# pick the priors for this set (gives NULL if no prior has been defined)
		set = sets[[i]]
		alpha = priors[[set]]$alpha 
		beta  = priors[[set]]$beta

		# solve 

		res <- RPA.iteration(S, sigma2.guess, epsilon, alpha, beta, sigma2.method, d.method)
		
		# Store the results	
		d.results[i,]<-res$d
		sigma2.results[[i]]<-res$sigma2
	
	}


	rpa.res = new("rpa", list(d=d.results, sigma2 = sigma2.results, cind=cind, sets = sets))

	# return result class
	rpa.res
}

