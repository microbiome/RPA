set.seed(4543)

require(RPA2)
require(RPA)
#fs <- list.files("~/Rpackages/RPA/github/RPA/pkg/R/", full.names = TRUE, pattern = ".R$"); for (f in fs) {source(f)}

set.seed(3453)
cels <- list.celfiles("CEL", full.names = T)[1:20]

eset.rpa1 <- RPA::rpa(cel.files = cels)
eset.rpa1b <- RPA::rpa.online(cel.files = cels)
eset.rpa2b <- RPA2::rpa.online(cel.files = cels)
eset.rpa2 <- RPA2::rpa(cel.files = cels)

e1 <- exprs(eset.rpa1)
e2 <- exprs(eset.rpa2)
e3 <- exprs(eset.rpa1b$expressionSet)
e4 <- exprs(eset.rpa2b$expressionSet)

par(mfrow = c(2,2))
rs <- sample(prod(dim(e1)), 1e4)
plot(as.vector(e1)[rs], as.vector(e2)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e2)[rs]))
plot(as.vector(e1)[rs], as.vector(e3)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e3)[rs]))
plot(as.vector(e1)[rs], as.vector(e4)[rs], pch = ".", main = cor(as.vector(e1)[rs], as.vector(e4)[rs]))

