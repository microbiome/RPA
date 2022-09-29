# Run 16.2.2013 with RPA_1.15.10

# ------------------------------------------------------------------

# Preprocess the data with Online-RPA 
require(RPA)
cels <- list.celfiles("CEL", full.names = T)
res <- RPA::rpa.online(cel.files = cels, batch.size = 200,                 		 
            save.batches.dir = ".", unique.run.identifier = "RPA-20130216", speedup = FALSE, mc.cores = 4)

# ------------------------------------------------------------------

# Investigate output

# Pick the applied parameters 
params <- res$params

# Pick the preprocessed expression matrix and customize sample names
emat <- exprs(res$expressionSet) 
colnames(emat) <- sapply(strsplit(colnames(emat), "\\."), function (x) {x[[1]]})
colnames(emat) <- sapply(strsplit(colnames(emat), "/"), function (x) {x[[2]]})

# Pick Hyperparameters
probe.parameters <- res$probe.parameters
probe.parameters.evolution <- res$probe.parameters.evolution

# Save the results
print("Saving all data")
save(res, emat, params, probe.parameters, probe.parameters.evolution, file = "LukkOutput201302016.RData", compress = "xz")

print("Saving probe parameters")
hgu133a.rpa.priors <- probe.parameters
save(hgu133a.rpa.priors, file = "HGU133A-RPA-priors.RData", compress = "xz")

# ----------------------------------------------------------------------



