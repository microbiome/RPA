run.id <- "RPA-run-id-0.821397874892181"

batch.size <- 200
datapath <- "~/Rpackages/RPA/OnlineLearning/"

# Hyperparameter files
fs <- list.files(datapath, pattern = run.id, full.names = TRUE)
fs <- fs[grep("hyper.RData", fs)]
o <- order(as.numeric(sapply(strsplit(fs, "-"), function (x) {x[[6]]})))
fs <- fs[o]

# Final estimated hyperparameters
final.hp <- paste(datapath, run.id, "-RPA-hyperparameters.RData", sep = "") # hyper.parameters
load(final.hp)

# Hyperparameter evolution
alphas <- c()
betamat <- NULL
for (f in fs) {
  load(f)
  alphas[[f]] <- alpha  
  betamat <- cbind(betamat, unlist(betas))
}
varmat <- t(t(betamat)/alphas)

# cumulative sample size from batches
cumulative.N <- batch.size * (1:length(fs))
cumulative.N[[length(cumulative.N)]] <- 5372 # 2 * alphas[[length(alphas)]] # alpha = N/2 + prior

# Plot parameters
#plot(cumulative.N, alphas, main = "alpha")
#hist(varmat[sample(length(varmat), 1e4)],100) # variance histogram
#hist(sqrt(varmat[sample(length(varmat), 1e4)]),100) # std histogram

# Final estimated variances
#final.vars <- unlist(hyper.parameters$variances)
#hist(sqrt(final.vars[sample(length(final.vars), 1e4)]), 1e2)

# Affinities
load("~/data/RPA/emat.rpa.online.storage.RData") # emat.rpa.online.storage
load("~/data/RPA/emat.rpa.online.storage.ENSG.RData") # emat.rpa.online.storage.ENSG

# emat.rpa.online.storage


# Probe std evolution
#probe <- sample(rownames(varmat), 1); plot(cumulative.N, sqrt(varmat[probe, ]), type = "b", xlab = "Cumulative sample size", ylab = "Standard deviation", main = probe, las = 1)


#############################################################

library(reshape)
library(ggplot2)
set <- "219877_at" # sample(rownames(varmat), 1)
#set <- "218054_s_at"
#mat <- sqrt(varmat[grep("^1007_s_at", rownames(betamat)),])
mat <- sqrt(varmat[grep(paste("^", set, sep = ""), rownames(betamat)),])
df <- melt(mat)
df$X2 <- cumulative.N[df$X2]
#ggplot(df) + aes(x = X2, y = value, group = X1) + geom_line()
p <- ggplot(df) + aes(x = X2, y = value, group = X1, col = X1) + geom_line()
jpeg("probestd.evolution.example.probeset.jpg")
print(p)
dev.off()