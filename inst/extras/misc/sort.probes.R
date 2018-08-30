# (C) Leo Lahti 2012. FreeBSD license.

# Install libraries
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(c("affy", "affydata", "RPA"))

# Load libraries
require(affy)
require(affydata)
require(RPA)

# Load example data
data(Dilution)

# Define the probesets to check
sets <- geneNames(Dilution)[1:2]

# Robust Probabilistic Averaging model
rpa.results <- RPA.pointestimate(Dilution, sets)

# Probe affinity effects
af <- unlist(lapply(sets, function (set) {rpa.results[[set]]$affinity}))

# Probe-specific noise (variance)
s2 <- unlist(lapply(sets, function (set) {rpa.results[[set]]$sigma2}))

# PM probe indices
pmind <- unlist(pmindex(Dilution)[sets])

# Probe effect table
df <- data.frame(list(pmindex = pmind, affinity = af, variance = s2))

# Probe effect table ordered by absolute affinity effects
df <- df[order(abs(df$affinity), decreasing = TRUE),]

# Investigate the output
print(head(df))
