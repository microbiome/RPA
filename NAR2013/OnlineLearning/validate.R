# ---------------------------------------------------------------------------

res2 <- rpa(cel.files = cels)
emat2 <- exprs(res2)[rownames(emat), ]
colnames(emat2) <- sapply(strsplit(colnames(emat2), "\\."), function (x) {x[[1]]})
emat2 <- emat2[rownames(emat), colnames(emat)]

load("~/data/RPA/Lukk-RPA.RData")
emat3 <- emat.rpa[rownames(emat),  colnames(emat)]

library(affy)
res4 <- rma(ReadAffy(filenames = cels))
emat4 <- exprs(res4)
colnames(emat4) <- sapply(strsplit(colnames(emat4), "\\."), function (x) {x[[1]]})
emat4 <- emat4[rownames(emat),  colnames(emat)]

###########################################################################

rs <- sample(prod(dim(emat)), 1e3)
pairs(cbind(rpaonline.new = emat[rs], rpa = emat2[rs], rpaonline.paper = emat3[rs], RMA = emat4[rs]))

