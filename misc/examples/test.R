library(RPA)

## Load example data set
require(ALLMLL)
data(MLL.A)

## Compute RPA for specific probesets
eset.rma <- rma(MLL.A)
#ps <- apply(exprs(eset.rma), 1, var)
#rev(sort(abs(ps)))[][1:10]

#fs <- list.files("~/local/Rpackages/RPA/SVN/RPA/R", full.names = TRUE, pattern = ".R$")
#for (f in fs) {source(f)}

set <- "214218_s_at";
rpa.estimate <- RPA.pointestimate(MLL.A, set, affinity.method = "rpa")
#cor(as.vector(rpa.estimate$d), exprs(eset.rma)[set,])
#pdf("~/tmp/tmp.pdf");rpa.plot(set, rpa.estimate);dev.off()
pdf("~/tmp/tmp.pdf");
rpa.plot(set, rpa.estimate, external.signal = exprs(eset.rma)[set,])
dev.off()


