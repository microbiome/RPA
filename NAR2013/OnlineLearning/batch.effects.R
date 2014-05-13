library(RPA)
library(ber)

# Load RPA-online data and class information vector
load("data/classInfo.RData") # classInfo
load("data/emat.rpa.online.storage.RData") # emat.rpa.online.storage

class.minsize <- 2
tab <- table(classInfo)
selected.classes <- names(which(tab >= class.minsize) )
s <- names(classInfo)[which(classInfo %in% selected.classes)]
cl <- classInfo[s]

# Sample batch vector
batch.vector <- factor(cl)
emat <- emat.rpa.online.storage[, names(cl)]
rm(emat.rpa.online.storage)

#########################################

remove.batches <- function (emat, cl) {

  # Remove batch effects with a linear model
  emat.adj <- matrix(NA, nrow = nrow(emat), ncol = ncol(emat))
  dimnames(emat.adj) <- dimnames(emat)
  for (i in 1:nrow(emat)) {
    print(i/nrow(emat))
    mod <- lm(emat[i,]~factor(cl))
    emat.adj[i,] <- summary(mod)$residuals + summary(mod)$coefficients["(Intercept)", "Estimate"]
  }

  emat.adj
} 


print("lm")
emat.rpa.online.lm <- remove.batches(emat, cl)
save(emat.rpa.online.lm, file = "data/emat.rpa.online.lm.RData")
rm(emat.rpa.online.lm)
gc()

# intercept + residuals = data - estimated.batch.term
#> summary(mod)$residuals[n] + summary(mod)$coefficients["(Intercept)", "Estimate"]
#GSM161547 GSM161548 GSM161549 
# 9.118724  9.145382  9.874999 
#> v[n] - (-0.2357360)
#GSM161547 GSM161548 GSM161549 
# 9.118724  9.145382  9.874999 

########################################


# Add batch-corrected versions of the data

print("bg")
emat.rpa.online.berbg   <- t(ber_bg(t(emat), batch.vector, partial = FALSE, nSim = 1e3))
save(emat.rpa.online.berbg, file = "data/emat.rpa.online.berbg.RData")
rm(emat.rpa.online.berbg)
gc()

print("mean")
emat.rpa.online.bermean <- t(mean_centering(t(emat), batch.vector))
save(emat.rpa.online.bermean, file = "data/emat.rpa.online.bermean.RData")
rm(emat.rpa.online.bermean)
gc()

print("sd")
emat.rpa.online.bersd <- t(standardization(t(emat), batch.vector))
save(emat.rpa.online.bersd, file = "data/emat.rpa.online.bersd.RData")
rm(emat.rpa.online.bersd)
gc()


# ------------------------------------------------------------------------

load("data/emat.rpa.online.storage.ENSG.RData") # 
emat <- emat.rpa.online.storage.ENSG[, names(cl)]
rm(emat.rpa.online.storage.ENSG)

print("lm")
emat.rpa.online.lm.ENSG <- remove.batches(emat, cl)
save(emat.rpa.online.lm.ENSG, file = "data/emat.rpa.online.lm.ENSG.RData")
rm(emat.rpa.online.lm.ENSG)
gc()

# Add batch-corrected versions of the data

print("bg")
emat.rpa.online.berbg.ENSG <- t(ber_bg(t(emat), batch.vector, partial = FALSE, nSim = 1e3))
save(emat.rpa.online.berbg.ENSG, file = "data/emat.rpa.online.berbg.ENSG.RData")
rm(emat.rpa.online.berbg.ENSG)
gc()

#print("mean")
emat.rpa.online.bermean.ENSG <- t(mean_centering(t(emat), batch.vector))
save(emat.rpa.online.bermean.ENSG, file = "data/emat.rpa.online.bermean.ENSG.RD#ata")
rm(emat.rpa.online.bermean.ENSG)
gc()

print("sd")
emat.rpa.online.bersd.ENSG   <- t(standardization(t(emat), batch.vector))
save(emat.rpa.online.bersd.ENSG, file = "data/emat.rpa.online.bersd.ENSG.RData")
gc()
