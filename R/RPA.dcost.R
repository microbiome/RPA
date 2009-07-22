RPA.dcost <-
function (d,sigma2,S) {
	# S : 	observed probe-level signals, T x P i.e. arrays x probes
	# d : 	assumed "real" signal for which this is a cost function
	# sigma2 : probe-specific variances
	M = S-d
	-sum((1/(2*sigma2))*((colSums(M)^2)/(nrow(S)+1) - colSums(M^2)))
}

