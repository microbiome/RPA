# 14.2.2013 

require(RPA2)
#fs <- list.files("~/Rpackages/RPA/github/RPA/pkg/R/", full.names = TRUE, pattern = ".R$")
#for (f in fs) {source(f)}

cel.list <- list.celfiles("~/data/RPA/CELarchive", full.names = T)
ns <- rev(c(1000, 5000, 10000, 20000))
bs <- 250

runtimes <- list()
for (n in ns) {
  print(n)
  cels <- sample(cel.list, n)
  runtimes[[as.character(n)]] <- system.time(eset.rpa.online <- RPA2::rpa.online(cel.files = cels, batch.size = bs, mc.cores = 6) )
  gc()
}

save(runtimes, ns, file = "runtimes-20130214.RData")

#source("runtimes-plot.R")

# mc.cores: 1 / batch.size: 50
#> runtimes
#$`300`
#    user   system  elapsed 
#2672.119   17.721 2712.142 

# mc.cores: 4 -> 2.3-fold speedup / batch.size: 50
#> runtimes
#$`300`
#    user   system  elapsed 
#2887.931   47.074 1173.002 

# mc.cores: 4 -> 1.56-fold speedup / batch.size: 150
#> runtimes
#$`300`
#    user   system  elapsed 
#1542.016   33.477  750.546 

# mc.cores: 4 / batch.size: 300
#$`300`
#    user   system  elapsed 
#1169.950   38.510  721.013 
# Same with ENSG:
#    user   system  elapsed 
#1543.276   31.689  745.067 


# batch.size `500`, mc.cores 4
#    user   system  elapsed 
# 1516.946  134.349 6657.467 


# Conclusion: 
# batch.size has remarkable but still smallish effect. 
# Larger batches are faster when they do not exhaust the memory.
# 4-times more cores will yield >2-fold speedup
# -> Use many cores and as large batches as possible so that the memory is not exhausted