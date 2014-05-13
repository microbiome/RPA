set.seed(4543)

require(RPA2)
require(RPA)
#fs <- list.files("~/Rpackages/RPA/github/RPA/pkg/R/", full.names = TRUE, pattern = ".R$"); for (f in fs) {source(f)}

cel.list <- sample(list.celfiles("CEL", full.names = T), 200)

ns <- 200
cels <- sample(cel.list, ns)
bs <- 100

eset.rpa.online.speedup <- RPA2::rpa.online(cel.files = cels, batch.size = bs, mc.cores = 4, speedup = TRUE, unique.run.identifier = "SPEEDUP", keep.batch.files = TRUE)
eset.rpa.online.nospeedup <- RPA2::rpa.online(cel.files = cels, batch.size = bs, mc.cores = 4, speedup = FALSE, unique.run.identifier = "NOSPEED", keep.batch.files = TRUE)
eset.rpa.online.nospeedup2 <- RPA2::rpa.online(cel.files = cels, batch.size = bs, mc.cores = 4, speedup = FALSE, unique.run.identifier = "NOSPEED2", keep.batch.files = TRUE, rseed = 13442)
eset.rpa.online.standard <- RPA2::rpa(cel.files = cels)
eset.rpa.online.standard.aff <- RPA2::rpa(cel.files = cels, summarize.with.affinities = TRUE)
eset.rpa.online.speedup.aff <- RPA2::rpa.online(cel.files = cels, batch.size = bs, mc.cores = 4, speedup = TRUE, unique.run.identifier = "SPEEDUP", keep.batch.files = TRUE, summarize.with.affinities = TRUE)

abatch <- ReadAffy(filenames = cels); eset2 <- RPA::rpa(abatch); e9 <- exprs(eset2); plot(as.vector(e1)[rs], as.vector(e9)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e9)[rs]))

eset.rpa.online.release <- RPA::rpa.online(cel.files = cels, batch.size = bs, mc.cores = 4, unique.run.identifier = "RELEASE", keep.batch.files = TRUE)

e1 <- exprs(eset.rpa.online.speedup$expressionSet)
e2 <- exprs(eset.rpa.online.nospeedup$expressionSet)
e3 <- exprs(eset.rpa.online.nospeedup2$expressionSet)
e4 <- exprs(eset.rpa.online.release$expressionSet)
e5 <- exprs(eset.rpa.online.standard)
e6 <- exprs(eset.rpa.online.standard.aff)
e7 <- exprs(eset.rpa.online.speedup.aff$expressionSet)


rs <- sample(prod(dim(e1)), 1e4)
par(mfrow = c(3,3))
plot(as.vector(e1)[rs], as.vector(e2)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e2)[rs]))
plot(as.vector(e1)[rs], as.vector(e3)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e3)[rs]))
plot(as.vector(e1)[rs], as.vector(e4)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e4)[rs]))
plot(as.vector(e1)[rs], as.vector(e5)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e5)[rs]))
plot(as.vector(e1)[rs], as.vector(e6)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e6)[rs]))
plot(as.vector(e1)[rs], as.vector(e7)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e7)[rs]))
plot(as.vector(e1)[rs], as.vector(e8)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e8)[rs]))
plot(as.vector(e1)[rs], as.vector(e10)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e10)[rs]))


