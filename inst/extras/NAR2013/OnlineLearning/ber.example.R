library(golubEsets)
library(vsn)
data(Golub_Merge)
E <- exprs(vsn2(Golub_Merge))
batch <- Golub_Merge$Source
BMPB <- Golub_Merge$BM.PB
BMPB <- data.frame(BMPB)
Eadj <- ber(t(E),batch,BMPB)
