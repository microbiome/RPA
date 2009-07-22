RPA.sigma2.update <-
function (d, S=S, alphahat, beta, sigma2.method="var") {
	# S : 	observed probe-level signals, T x P i.e. chips x probes
	# d : estimated probe-set level differential expression profile
	# alphahat, beta : parameters for inverse Gamma distribution
 
	# arrays x probes matrix; observations vs. estimated real signal
	R<-S-d 

	if (sigma2.method == "var") {
	        # Assume uninformative priors alpha, beta -> 0
		# NOTE: our formulation converges to variance with large sample sizes
		s2 = apply(R,2,function(r){var(r)}) 
	} else if (sigma2.method == "basic") {
		# update betahat
		betahat <- .5*sapply(seq(ncol(R)),function(i){r=R[,i];(2*beta[[i]] + sum(r^2) - sum(r)^2 / (nrow(S)+1))})

		# updated sigma2
		s2 = betahat / (alphahat + 1) # = sigma2 = mode for invgam(sig^2 | alphahat,betahat)
	}
	# return updated sigma2
	s2
}

