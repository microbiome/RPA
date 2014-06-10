library(RPA)
#fs <- list.files("../pkg/R", full.names = TRUE, pattern = ".R$")
#for (f in fs) {source(f)}

cels <- sample(list.files("CEL", full.names = TRUE), 300)

res1 <- rpa(cel.files = cels) 
#res2 <- rpa.online(cel.files = cels, batch.size = Inf, save.batches = TRUE)
#res3 <- rpa.online(cel.files = cels, batch.size = Inf, save.batches = FALSE)
#res4 <- rpa.online(cel.files = cels, batch.size = 20, save.batches = FALSE)
res5 <- rpa.online(cel.files = cels, batch.size = 100, save.batches = TRUE, save.batches.dir = "~/tmp/", mc.cores = 4)



emat1 <- exprs(res1)
#emat2 <- exprs(res2$expressionSet)
#emat3 <- exprs(res3$expressionSet)
#emat4 <- exprs(res4$expressionSet)
emat5 <- exprs(res5$expressionSet)

rs <- sample(prod(dim(emat1)), 1e4)
#pairs(cbind(emat1[rs], emat2[rs], emat3[rs], emat4[rs], emat5[rs]))
pairs(cbind(emat1[rs], emat5[rs]))

#save.image()