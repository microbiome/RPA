# Test vignette with function updates before build & check
# ~/bin/R-2.13.0/bin/R
library(RPA)
fs <- list.files("../RPA/R/", full.names = TRUE)
for (f in fs) {source(f)}
Sweave("RPA.Rnw")
