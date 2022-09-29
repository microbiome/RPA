# Run 12.11.2012 with RPA_1.15.10

# ------------------------------------------------------------------

# Preprocess the data with Online-RPA 
require(RPA)
res <- rpa.online(cel.path = "CEL", batch.size = 200,                 		 
            save.batches.dir = ".", unique.run.identifier = "RPA-20121104")

# ------------------------------------------------------------------

# Investigate output

# Pick the applied parameters 
params <- res$params

# Pick the preprocessed expression matrix and customize sample names
emat <- exprs(res$expressionSet) 
colnames(emat) <- sapply(strsplit(colnames(emat), "\\."), function (x) {x[[1]]})
colnames(emat) <- sapply(strsplit(colnames(emat), "/"), function (x) {x[[2]]})

# Pick Hyperparameters
probe.parameters <- res$hyper.parameters
probe.parameters2 <- probetable(res$hyper.parameters)
hyper.parameters.evolution <- res$hyper.parameters.evolution

# -----------------------------------------------------

# Save the results

save(res, emat, params, probe.parameters, hyper.parameters.evolution, file = "LukkOutput20121104.RData", compress = "xz")


# load("LukkOutput20121104.RData")
library(RPA);
probe.parameters <- probetable(probe.parameters)
save(emat, probe.parameters, file = "Online-RPA-Lukk2010-Affy.RData", compress = "xz")

# Expression matrix for publication
write.table(round(emat, 3), file = "HGU133A-Lukk2010-Affy-RPA-ExpressionMatrix.tab", sep = "\t", quote = FALSE, row.names = FALSE)


# Hyperparameter table for publication
probe.parameters.rounded <- probe.parameters
probe.parameters.rounded[, c("alpha", "beta", "tau2", "mu")] <- round(probe.parameters[, c("alpha", "beta", "tau2", "mu")], 3)
write.table(probe.parameters.rounded, file = "HGU133A-Lukk2010-Affy-RPA-ProbeParams.tab", sep = "\t", quote = FALSE, row.names = FALSE)



#load("~/data/RPA/Lukk-RPA.RData")
#emat3 <- emat.rpa[rownames(emat),  colnames(emat)]
#rs <- sample(prod(dim(emat)), 1e5)
#plot(emat[rs], emat3[rs], pch = '.')



