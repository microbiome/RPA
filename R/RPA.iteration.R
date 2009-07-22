RPA.iteration <-
function(S, sigma2.guess, epsilon = 10^(-3), alpha=NULL, beta=NULL, sigma2.method = "var", d.method = "fast") {
	# sigma2.guess : initial guess for probe-specific variances
	# epsilon : convergence threshold
	# S : data; arrays x probes matrix of probe-wise fold-changes
	# alpha,beta : priors for the inverse Gamma distribution
	# 	       used only when sigma2.method = "basic"
	#	       These are vectors (each element corresponds to one probe)
	# sigma2.method: optimization method for sigma2
	# 		 "basic": optimization with any priors
	# 		 "var":   uses the fact that the cost function converges
	#		 	  to variance with larger sample sizes
	# d.method: "fast": weighted mean 
	# 	    	    over the probes, weighted by probe variances
	#		    The solution converges to this with large sample size
	# 	     "basic": optimization scheme to find a mode
	#	     	      used in the original publication; slow


	# if no prior has been given, use noninformative default priors
	P = length(sigma2.guess) # number of probes
	if (length(alpha)==0) {alpha = rep(1e-6,P)}
	if (length(beta)==0) {beta = rep(1e-6,P)}

	T <- nrow(S) # Number of arrays (except control)
	alphahat <- T/2 + alpha

	if (d.method == "fast") {converged = TRUE} else {converged = FALSE}

	#print(" Alternate optimization of d and sigma2 until convergence")
	sigma2 <- sigma2.guess #This initial value is actually used in the iterations
	#This initial value is not used in computation, just for checking convergence at first iteration
	sigma2.old <- sigma2.guess + 1000*epsilon 

	d<-rep(max(S),T) #This initial value is not used in computation, just for checking convergence at first iteration
	d.old<-(-d) #This initial value is not used in computation, just for checking convergence at first iteration
	d.init<-rep(0,T) #initial value for d in the optimization 

	cnt<-0

	# optimize until convergence
	loopcnt = 0
	while ((max(abs(c(sigma2-sigma2.old,d-d.old)))>epsilon | !converged) && loopcnt < 1e6) {
		cnt<-cnt+1	
		#print(cnt)	

		d.old<-d
		sigma2.old<-sigma2

		# optimize d
		# start optimization from zero signal at each iteration
		# (alternatively, from d.old i.e. the previous estimate)
		if (d.method == "basic") {
		  res <- optim(d.init,fn=RPA.dcost,method="BFGS",sigma2=sigma2,S=S)

		  # Check convergence at this round
		  converged = (res$convergence==0)

    		  if (converged) {
			# update d
			d<-res$par
		  } else {
			# following iteration count to avoid 
			# 	potentially infinite loops
		    	loopcnt = loopcnt + 1 

			# If convergence problems occur, initialize with the mean signal of the probes
			# also do some random perturation in case we would end up here twice
			# This is very rare according to our experience, occuring less than in 1 / 50000 
			# analyzed probesets 
			d.init<-rowMeans(S)
			d.init <- rnorm(length(d.init),mean = d.init, sd = sd(d.init))
			print("Initialized d with perturbed probeset mean")

		}
	       } else if (d.method == "fast") {
	       	 d <- d.update.fast(S, sigma2)
	       }

		#update sigma2
		sigma2 <- RPA.sigma2.update(d,S=S,alphahat=alphahat,beta=beta, sigma2.method)
	}

	list(d=d,sigma2=sigma2)
}

