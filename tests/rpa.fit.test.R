library(RPA)

# Generating toy data
P <- 10 # number of probes (observations)
n <- 20 # number of samples in the multivariate observation
probe.noise <- rchisq(P, 1)
probe.affinity <- rnorm(P, sd = probe.noise)
mu.real <- 2
signal <- rnorm(n)
dat <- array(NA, dim = c(P, n))
for (p in 1:P) {dat[p, ] <- rnorm(length(signal), mean = signal + mu.real + probe.affinity[[p]], sd = sqrt(probe.noise[[p]]))}

# Fit the RPA model
res <- rpa.fit(dat)

#pdf("~/tmp/tmp.pdf")
#rpa.plot(dat, rpa.fit.object = res, plots = "data")
#dev.off()

# Compare estimates to truth
#pdf("~/tmp/tmp.pdf")
par(mfrow = c(2,2))
rpa.plot(dat, rpa.fit.object = res, plots = "data", main = "Data and summary")
barplot(rbind(signal + mu.real, res$mu), beside = TRUE, main = "Signal", xlab = "Sample", ylab = "Signal")
barplot(rbind(probe.affinity, res$affinity), beside = TRUE, main = "Affinity", xlab = "Probe index", ylab = "Affinity")
barplot(rbind(probe.noise, res$sigma2), beside = TRUE, main = "Probe variance", xlab = "Probe index", ylab = "Variance")
dev.off()

#rpa.plot(dat, rpa.fit.object = res)
