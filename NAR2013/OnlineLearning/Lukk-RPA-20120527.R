# Tested 27.5.2012. OK.

require(RPA)
fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}

cels <- sample(list.celfiles("CEL", full.names = T))

runtime <- system.time(emat.rpa.online <- exprs(rpa.online(cel.files = cels, priors = list(alpha = 1, beta = 1), batch.size = 200, save.batches = TRUE, save.batches.dir = "Lukk-RPA-20120527") ))

colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "/"), function (x) {x[[2]]})
colnames(emat.rpa.online) <- sapply(strsplit(colnames(emat.rpa.online), "\\."), function (x) {x[[1]]})

# Compare to earlier:
load("~/data/RPA/Lukk-RPA.RData")
emat.rpa.online <- emat.rpa.online[, colnames(emat.rpa)]
#inds <- sample(prod(dim(emat.rpa)),1e4); cor(emat.rpa[inds], emat.rpa.online[inds])
# Confusing -> rather poor correlation 0.96??

save(runtime, emat.rpa.online, file = "emat.rpa.online.RData", compress = "xz")

#####################################

# Tested 27.5.2012. OK.
require(RPA)
fs <- list.files("~/Rpackages/RPA/RPA/R/", full.names = TRUE, pattern = ".R$")
for (f in fs) {source(f)}
cels <- sample(list.celfiles("CEL", full.names = T))
emat.rpa.online.test <- exprs(rpa(cel.files = sample(cels, 300)))
colnames(emat.rpa.online.test) <- sapply(strsplit(colnames(emat.rpa.online.test), "/"), function (x) {x[[2]]})
colnames(emat.rpa.online.test) <- sapply(strsplit(colnames(emat.rpa.online.test), "\\."), function (x) {x[[1]]})
save(emat.rpa.online.test, file = "emat.rpa.online.test.RData", compress = "xz")

load("emat.rpa.online.RData")

rind <- sample(prod(dim(emat.rpa.online.test)), 1e4)
s <- colnames(emat.rpa.online.test)
plot(emat.rpa.online[, s][rind], emat.rpa.online.test[, s][rind])
cor(emat.rpa.online[, s][rind], emat.rpa.online.test[, s][rind])
