#~/bin/R-2.13.0/bin/R
require(affy)
require(MCMCpack)

require(RPA)
set.seed(11122)
# Generate and fit toydata, learn hyperparameters
fs <- list.files("~/Rpackages/RPA/RPA/R", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}

# Generate random probeset data
P <- 11   # number of probes
N <- 2400 # number of arrays
real <- sample.probeset(P = P, n = N, shape = 3, scale = 1, mu.real = 4)
dat <- real$dat # probes x samples

#######################################################

start <- Sys.time()
print(system.time(source("test3.R")))
et <- Sys.time() - start
print(et)

# Check:
#par(mfrow = c(2,1))
# Online
plot(sqrt(real$variance), sqrt(s2), main = cor(sqrt(real$variance), sqrt(s2))); abline(0,1)
# Total
#plot(sqrt(real$variance), sqrt(s2.tot), main = cor(sqrt(real$variance), sqrt(s2.tot))); abline(0,1)
# Excellent; moreover, online version beats complete data version


# TODO: speedups


