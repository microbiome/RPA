#~/bin/R-2.13.0/bin/R
require(affy)
require(RPA)

set.seed(11122)
# fs <- list.files("../RPA/R", full.names = TRUE, pattern = ".R$")
# for (f in fs) {source(f)}

# Generate and fit toydata, learn hyperparameters
P <- 11   # number of probes
N <- 20 # number of arrays
real <- sample.probeset(P = P, n = N, shape = 3, scale = 1, mu.real = 4)
dat <- real$dat # probes x samples

#######################################################

start <- Sys.time()

alpha0 <- beta0 <- 1e-2

# Set priors
alpha <- alpha0
beta  <- rep(beta0, P)
s2 <- beta/alpha

# Operate in batches

step <- 10
for (ni in seq(1, N, step)) {
  batch <- ni:(ni+step-1)  
  # Estimate s2 with EM-type procedure
  my.dat <- dat[,batch]
  s2 <- s2.update(my.dat, alpha, beta, s2.init = beta/alpha, th = 1e-2)

  #beta.test <- update.beta.EM(my.dat, alpha, beta, th = 1e-2)

  # Update hyperparams: 
  alpha <- update.alpha(ncol(my.dat), alpha)
  beta <- alpha * s2 # from the mode of s2, i.e. beta/alpha

}

# Final variance
s2 <- beta/alpha

et <- Sys.time() - start
print(et)

# Online
par(mfrow=c(2,1))
plot(sqrt(real$variance), sqrt(s2), main = cor(sqrt(real$variance), sqrt(s2))); abline(0,1)
plot(real$mu.real + real$d, d.update.fast(dat, s2)); abline(0,1)


