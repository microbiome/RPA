d.update.fast <-
function (S, sigma2) {
	 # With large sample sizes when T->Inf
	 # the d converges to the weighted mean 
	 # over the probes, weighted by	probe variances	 
	 rowSums(t(  (t(S)/sigma2) / (sum(1/sigma2))))
}

