# Run 28.5.2012 with RPA_1.13.02
#hgu133ahsensgcdf_15.1.0
# 6.6. with hgu133ahsensgcdf_14.1.0 for compatibility

require(RPA)
#fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names = TRUE, pattern = ".R$")
#for (f in fs) {source(f)}

set.seed(23)
cels <- sample(list.celfiles("CEL", full.names = T))

eset <- rpa.online(cel.files = cels, priors = list(alpha = 1, beta = 1), batch.size = 200,                 		 
            save.batches = TRUE, save.batches.dir = ".", keep.batch.files = TRUE, cdf = "HGU133A_Hs_ENSG")

	    	   
emat.rpa.online <- exprs(eset) 
colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "\\."), function (x) {x[[1]]})
colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "/"), function (x) {x[[2]]})

emat.rpa.online.storage.ENSG <- emat.rpa.online 
save(emat.rpa.online.storage.ENSG, file = "emat.rpa.online.storage.ENSG.RData", compress = "xz")


#print("RPA4") # original
#sets <- rownames(emat.rpa.online.storage)
#s <- colnames(emat.rpa.online.storage)
#load("~/data/RPA/Lukk-RPA.RData")
#emat.rpa.online.usedinexperiments <- emat.rpa[sets, s]

#rsample <- sample(prod(dim(emat.rpa.online.usedinexperiments)), 1e4)
#s <- colnames(emat.rpa.online.usedinexperiments)
#tab <- cbind(RPA.online.init = emat.rpa.online.usedinexperiments[sets,s][rsample], RPA.storage = emat.rpa.online.storage[sets,s][rsample])
#print(cor(tab, use = "pairwise.complete.obs"))
#pairs(tab)



