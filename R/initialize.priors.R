initialize.priors <-
function (abatch, sets, alpha = 1e-6, beta = 1e-6, d = NULL) {

	#setClass("rpa.priors", contains = "list", where=topenv(parent.frame()))

	# set priors for d and sigma2
	# skip in the current version the prior for d

	priors = vector(length=length(sets),mode="list")

	for (set in sets) {	
		nprobes = nrow(pm(abatch,set))
		alphas = rep(alpha,nprobes)
		betas = rep(beta,nprobes)
		ds = d 
		priors[[set]] = list(alpha = alphas, beta=betas, d = ds) 
	}

	new("rpa.priors", priors)
}

