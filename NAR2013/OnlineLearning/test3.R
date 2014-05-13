alpha0 <- beta0 <- 1e-2

# Set priors
alpha <- alpha0
beta  <- rep(beta0, P)
s2 <- beta/alpha

# Operate in batches

step <- 200
for (ni in seq(1, N, step)) {

  print(ni)

  batch <- ni:(ni+step-1)
  datc <- t(centerData(t(dat[,batch])))

  # Estimate s2 with EM-type procedure using MCMC sampling for s2
  #s2 <- s2.EM(datc, alpha, beta, s2.init = s2, th = 1e-2, optim.method = "BFGS")
  s2 <- s2.update(datc, alpha, beta, s2.init = s2, th = 1e-2)

  # Update hyperparams: 
  alpha <- update.alpha(ncol(datc), alpha)
  beta <- alpha * s2 # from the mode of s2, i.e. beta/alpha

}

# Final variance
s2 <- beta/alpha

# Compare to that from whole data
#datc <- t(centerData(t(dat)))
#s2.tot <- s2.EM(datc, alpha0, rep(beta0, P), th = 1e-3, optim.method = "BFGS")


