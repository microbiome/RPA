# LL 7.2.2013
# Compare SVN release and experimental github versions
# See also RPA/OnlineLearning/comparison.R

library(RPA)
library(RPA2)

# Specify CEL files
cels <- list.celfiles("~/data/Kettunen13-LungCELs/Affy-data", full.names = T)
cels <- cels[-grep("HWik_z23T_110803.CEL", cels)]
cels <- cels[-grep("Hwik_z252_Barray_scan2_041203.CEL", cels)]
cels <- cels[-grep("HWik_z45N_300603.CEL", cels)]
cels <- cels[-grep("HWik_z45T_300603.CEL", cels)]
cels <- cels[-grep("HWik_z80T_310703.CEL", cels)]
cels <- cels[-grep("HWik_z80T_31scan2.CEL", cels)] 
cels <- cels[-grep("z170\\(2\\)-040203-HWik.CEL", cels)]


print("Release")
tt1 <- system.time(eset1 <- RPA::rpa(cel.files = cels))

print("RPA online singlebatch save")
tt3 <- system.time(eset3 <- RPA::rpa.online(cel.files = cels, batch.size = Inf, save.batches = TRUE, save.batches.dir = "."))

print("RPA online twobatch save")
tt7 <- system.time(eset7 <- RPA::rpa.online(cel.files = cels, batch.size = 20, save.batches = TRUE, save.batches.dir = "."))

# ---------------------------------------

print("Devel")
tt2 <- system.time(eset2 <- RPA2::rpa(cel.files = cels))

# ---------------------------------------

#fs <- list.files("RPA/pkg/R/", full.names=T); for (f in fs) {source(f)}

# RPA2 online -> almost 10-fold speedup !

print("RPA2 online singlebatch save")
tt4 <- system.time(eset4 <- RPA2::rpa.online(cel.files = cels, batch.size = Inf, save.batches.dir = "."))

# Check later
print("RPA2 online twobatch save")
tt8 <- system.time(eset8 <- RPA2::rpa.online(cel.files = cels, batch.size = 20, save.batches.dir = "."))

# -----------------------------------------

library(affy)
mat <- cbind(
    eset1 = as.vector(exprs(eset1)),
    eset2 = as.vector(exprs(eset2)),
    eset3 = as.vector(exprs(eset3$expressionSet)),
    eset7 = as.vector(exprs(eset7$expressionSet)),
    eset4 = as.vector(exprs(eset4$expressionSet)),
    eset8 = as.vector(exprs(eset8$expressionSet)) # check later
    )

print(cor(mat))

rs <- sample(prod(dim(eset1)), 1e4)
pairs(mat[rs,], pch = ".")

save.image()
 
tt1
tt2
tt3
tt4
tt7
tt8


# Also test affinity effect


#e4 <- apply(exprs(eset4$expressionSet), 1, function (x) {x - mean(x)})
#e6 <- apply(exprs(eset6$expressionSet), 1, function (x) {x - mean(x)})
#plot(as.vector(e4), as.vector(e6), pch = ".")
#plot(as.vector(exprs(eset4$expressionSet)), as.vector(exprs(eset6$expressionSet)), pch = ".")
