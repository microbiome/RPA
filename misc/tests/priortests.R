set.seed(4543)

require(RPA)
require(RPA2)
fs <- list.files("~/Rpackages/RPA/github/RPA/pkg/R/", full.names = TRUE, pattern = ".R$"); for (f in fs) {source(f)}
#cel.list <- list.celfiles("~/data/RPA/CELarchive", full.names = T)
source("../KettunenCEL.R"); cel.list <- cels

ns <- 40
bs <- 20
cels <- sample(cel.list, ns)

tt0 <- system.time(res0 <- RPA::rpa(cel.files = cels))
tt1 <- system.time(res1 <- rpa.complete(cel.files = cels))
tt2 <- system.time(res2 <- rpa(cel.files = cels, probe.parameters = res1$probe.parameters))
tt2b <- system.time(res2b <- rpa(cel.files = cels, probe.parameters = res1$probe.parameters))
tt2c <- system.time(res2c <- rpa(cel.files = cels))
tt3 <- system.time(res3 <- eset.rpa.online.speedup <- rpa.online(cel.files = cels, batch.size = bs, mc.cores = 4, speedup = TRUE, unique.run.identifier = "SPEEDUP"))
tt4 <- system.time(res4 <- eset.rpa.online.speedup <- rpa.online(cel.files = cels, batch.size = bs, mc.cores = 4, speedup = TRUE, unique.run.identifier = "SPEEDUP2", probe.parameters = res3$probe.parameters))
tt5 <- system.time(res5 <- eset.rpa.online.speedup <- rpa.online(cel.files = cels, batch.size = bs, mc.cores = 4, speedup = TRUE, unique.run.identifier = "SPEEDUP2", summarize.with.affinities = TRUE))

rs <- sample(prod(dim(exprs(res0))), 1e4)
par(mfrow = c(3, 3))
plot(exprs(res0)[rs], exprs(res1$eset)[rs], pch = ".")
plot(exprs(res0)[rs], exprs(res2$eset)[rs], pch = ".")
plot(exprs(res0)[rs], exprs(res3$expressionSet)[rs], pch = ".")
plot(exprs(res0)[rs], exprs(res4$expressionSet)[rs], pch = ".")
plot(exprs(res0)[rs], exprs(res5$expressionSet)[rs], pch = ".")


# range(abs(exprs(res2) - exprs(res2b))) # Should be 0
# range(abs(exprs(res2) - exprs(res2c))) # Should have some variation