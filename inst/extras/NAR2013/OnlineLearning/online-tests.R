# Calculates correlation of the estimated
# signal and real signal in online-learning and
# alternatives.
# When batch size grows in online-learning ("step")
# the result converges to RPA calculated from the complete data
# When batch size is small compared to the data, 
# online learning is not better than median

#~/bin/R-2.13.0/bin/R
require(affy)
library(affydata)
data(Dilution)
require(RPA)

set.seed(11122)

# Generate and fit toydata, learn hyperparameters
fs <- list.files("../RPA/R", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}



# Generate random probeset
P <- 11   # number of probes
N <- 1000 # number of arrays

real <- sample.probeset(P = P, n = N, shape = 1, scale = 1, mu.real = 3)

dat <- real$dat # probes x samples


#beta.EM(dat, s2.old = rep(1, nrow(dat)), alpha0 = 1e-3, beta0 = 1e-3)

#############################################

# Split in pieces
alpha <- alpha.prior <- rep(3, P)
beta <- beta.prior <- rep(1, P)
alphas <- betas <- vars <- NULL
step <- 1000 # 500: online ~ median; 2000: online ~ total
ns <- floor(c(1, seq(step,N,step)))
for (i in 2:length(ns)) {

  #print(ns[[i]])
    
  # Pick the next bunch of samples
  n.start <- ns[[i-1]]
  n.stop <- ns[[i]]
  samples <- (n.start+1):n.stop

  # Train with updated priors
  estimated <- rpa.fit(dat[, samples], alpha = alpha, beta = alpha)

  # Update alpha and beta
  alpha <- estimated$alpha
  beta <- estimated$beta

  # Store hypers
  alphas <- rbind(alphas, alpha)
  betas <- rbind(betas, beta)

  # Estimated variances
  vars <- rbind(vars, estimated$sigma2)

}

########################################################

# Pick online-learned priors
online.alpha <- alpha
online.beta <- beta
online.s2 <- beta/alpha
online.d <- d.update.fast(dat, online.s2)

#############################################

# Model for complete data
estimated <- rpa.fit(dat, alpha = 3, beta = 1, d.method = "fast", sigma2.method = "fast")

# Update alpha and beta from modeling of total data set
total.alpha <- estimated$alpha
total.beta <- estimated$beta
total.s2 <- estimated$sigma2
total.d <- estimated$mu
mean.d <-colMeans(dat)
median.d <-apply(dat,2,median)

#plot(sqrt(estimated$sigma2), sqrt(estimated$beta/estimated$alpha));abline(0,1)

#############################################

# Take weighted average of the probes, given the online variances
d.real <- real$mu.real + real$d

tab <- cbind(real = d.real, 
            online = online.d,
	    total = total.d,
	    mean = mean.d,
	    median = median.d)

print(sort(cor(tab)[1,]))


#pairs(cbind(real = d.real, 
#            online = online.d,#
#	    total = total.d,
#	    mean = mean.d
#))


#plot(total.d[1:10], mean.d[1:10]); abline(0,1)

############################################

#par(mfrow = c(2,2))
#plot(ns[-1], vars[,1], main = "var convergence"); abline(real$variance[[1]], 0)
#plot(ns[-1], betas[,1]/alphas[,1], main = "var.ab convergence"); abline(real$variance[[1]], 0)

#sim.online <- sum((online.d - d.real)^2)

# Std estimate with total data in very good agreement with ground truth
# Online version has some problem
#vars <- cbind(online = sqrt(online.beta/online.alpha),
#      real = sqrt(real$variance),
#      total = sqrt(total.beta/total.alpha))
#pairs(vars, ylim = range(vars), xlim = range(vars))

##############################################

