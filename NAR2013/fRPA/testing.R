library("aroma.affymetrix");
cdf <- AffymetrixCdfFile$byChipType("HG-U133_Plus_2");
#cs <- AffymetrixCelSet$byName("Affymetrix-HeartBrain", cdf=cdf);
#cs <- AffymetrixCelSet$byName("Affymetrix-HeartBrain", cdf=cdf);

cs <- "/home/lei/data/RPA/LukkAtlas/" #sample(list.celfiles("/home/lei/data/RPA/LukkAtlas/", full.names = T), 100)

# RMA background correction
bc <- RmaBackgroundCorrection(cs);
csB <- process(bc);

# RMA quantile normalization
qn <- QuantileNormalization(csB, typesToUpdate="pm");
csN <- process(qn);
