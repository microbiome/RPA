cels <- list.celfiles("~/tmp/CEL/CEL", full.names = T)[1:20]
abatch <- ReadAffy(filenames = cels)

# Retrieve probe positions
# The indices correspond to the output of pm()
pN <- probeNames(abatch)
set.inds <- split(1:length(pN), pN) # pNList


# Generate and fit toydata, learn hyperparameters
P <- length(set.inds[[1]])   # number of probes
N <- length(cels) # number of arrays
real <- sample.probeset(P = P, n = N, shape = 3, scale = 1, mu.real = 4)
dat <- real$dat # probes x samples

# Put known data in the first probeset
pm(abatch)[set.inds[[1]],] <- exp(dat)


