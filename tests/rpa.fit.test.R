library(RPA)
#fs <- list.files("../RPA/R", full.names = TRUE, pattern = ".R$")
#for (f in fs) {source(f)}

# Generate random probeset
P <- 10
n <- 100
real <- sample.probeset(P = P, n = n, shape = 3, scale = 1, mu.real = 3)
dat <- real$dat

# Fit the RPA model for the probeset, start with uninformative priors
# small priors ~ 0 are uninformative
alpha <- rep(1 + 1e-3, P)
beta <- rep(1e-3, P)
estimated <- rpa.fit(dat, alpha = alpha, beta = beta)

# Check the results
rpa.plot(dat, rpa.fit.object = estimated, toydata.object = real, plots = "toydata")

# FIXME: compare sigma2.method robust and fast performance numerically
