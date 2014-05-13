# Tested 25.5.2012. OK.

require(RPA)
fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}

set.seed(23)
cels <- sample(list.celfiles("CEL", full.names = T))#[1:300]

#print("RMA")
#emat.rma <- exprs(rma(ReadAffy(filenames = cels)))

#nsets <- 1000 # 
#sets <- sample(rownames(emat.rma), nsets)
#nsets <- nrow(emat.rma)
#sets <- rownames(emat.rma)

batch.size <- 200 

#print("RPA1") # with temporary data storage
#runtimes[["rpa-online"]] <- system.time(emat.rpa.online <- exprs(rpa.online(cel.files = cels, priors = list(alpha = 1, beta = 1), batch.size = batch.size, sets = sets) ))
#                cel.path = NULL;
#		cel.files = cels; 
#                    sets = sets;
#                     cdf = NULL; 
#               bg.method = "rma";                              
#                  priors = list(alpha = 1, beta = 1);
#                 epsilon = 1e-2;
#                mc.cores = 1;
#                 verbose = TRUE;                          
#                 shuffle = TRUE;                          
#		 batch.size = batch.size
#                 batches = NULL; 
#          quantile.basis = NULL; 
#            save.batches = TRUE;
#	    save.batches.dir = "."; 
#	    keep.batch.files = TRUE; 
#	    unique.run.identifier = NULL
#source("rpa.online.debug.R") 

eset <- rpa.online(cel.files = cels, 
                  priors = list(alpha = 1, beta = 1),
		 batch.size = batch.size,
            save.batches = TRUE,
	    save.batches.dir = ".", 
	    keep.batch.files = TRUE)

emat.rpa.online <- exprs(eset) 

colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "\\."), function (x) {x[[1]]})
colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "/"), function (x) {x[[2]]})

emat.rpa.online.storage <- emat.rpa.online 
save(emat.rpa.online.storage, file = "emat.rpa.online.storage.RData", compress = "xz")



skip <- TRUE
if (skip) {
  print("RPA2") # without temporary data storage
                cel.path = NULL;
		cel.files = cels; 
                    sets = sets;
                     cdf = NULL; 
               bg.method = "rma";                              
                  priors = list(alpha = 1, beta = 1);
                 epsilon = 1e-2;
                mc.cores = 1;
                 verbose = TRUE;                          
                 shuffle = TRUE;                          
		 batch.size = batch.size
                 batches = NULL; 
          quantile.basis = NULL; 
            save.batches = FALSE;
	    save.batches.dir = "."; 
	    keep.batch.files = FALSE; 
	    unique.run.identifier = NULL
source("rpa.online.debug.R") 
emat.rpa.online2 <- exprs(eset)

print("RPA3")
emat.rpa.standard <- exprs(rpa(cel.files = cels, priors = list(alpha = 1, beta = 1), verbose = TRUE, sets = sets))

}


print("RPA4") # original
sets <- rownames(emat.rpa.online.storage)
s <- colnames(emat.rpa.online.storage)
load("~/data/RPA/Lukk-RPA.RData")
emat.rpa.online.usedinexperiments <- emat.rpa[sets, s]


#colnames(emat.rpa.online.usedinexperiments) 
#colnames(emat.rpa.online.storage)
#colnames(emat.rma)


#colnames(emat.rma) <- sapply(strsplit(colnames(emat.rma), "\\."), function (x) {x[[1]]})
#colnames(emat.rpa.standard) <- sapply(strsplit(colnames(emat.rma), "\\."), function (x) {x[[1]]})

#colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "\\."), function (x) {x[[1]]})
#colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "/"), function (x) {x[[2]]})

#colnames(emat.rpa.online2) <- sapply(strsplit(colnames(emat.rpa.online2), "\\."), function (x) {x[[1]]})
#colnames(emat.rpa.online2) <- sapply(strsplit(colnames(emat.rpa.online2), "/"), function (x) {x[[2]]})


rsample <- sample(prod(dim(emat.rpa.online.usedinexperiments)), 1e4)
s <- colnames(emat.rpa.online.usedinexperiments)
#tab <- cbind(RMA = emat.rma[sets,s][rsample], RPA.online.init = emat.rpa.online.usedinexperiments[sets,s][rsample], RPA.online = emat.rpa.online2[sets,s][rsample], RPA.storage = emat.rpa.online[sets,s][rsample], RPA.standard = emat.rpa.standard[sets,s][rsample])
tab <- cbind(RPA.online.init = emat.rpa.online.usedinexperiments[sets,s][rsample], RPA.storage = emat.rpa.online.storage[sets,s][rsample])
print(cor(tab, use = "pairwise.complete.obs"))

pairs(tab)



