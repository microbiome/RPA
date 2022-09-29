
require(RPA)
#fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names = TRUE, pattern = ".R$")
#for (f in fs) {source(f)}

set.seed(23)
cels <- sample(list.celfiles("CEL", full.names = T))

eset <- rpa.online(cel.files = cels, priors = list(alpha = 1, beta =
1), batch.size = 20, save.batches = TRUE, save.batches.dir = ".",
keep.batch.files = TRUE, unique.run.id = "uniq")

	    	   
emat.rpa.online <- exprs(eset) 
colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "\\."), function (x) {x[[1]]})
colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "/"), function (x) {x[[2]]})

#save(emat.rpa.online, file = "emat.rpa.online.RData", compress = "xz")

#############################################################

# Robust Probabilistic Averaging model
abatch <- ReadAffy(filenames = cels)
rpa.results <- RPA.pointestimate(abatch)
# Provide a table of probe affinity and variance parameters
probe.parameter.table <- probe.performance(rpa.results)

##############################################################

# Load online estimates
load("uniqRPA-affinities.RData")
load("uniqRPA-hyperparameters.RData")

#########################################################

# Compare standard and online estimates

#hyper.parameters

df <- data.frame(pm.index = names(unlist(affinities)), affinity = unlist(affinities), variance = unlist(lapply(hyper.parameters$betas, function (beta) {beta/hyper.parameters$alpha})))

plot(probe.parameter.table$affinity, df$affinity); abline(0,1)

plot(probe.parameter.table$variance, df$variance); abline(0,1)