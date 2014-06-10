load("LukkOutput201302016.RData")
emat.new <- emat

load("~/data/RPA/Lukk-RPA.RData")
emat.old <- emat.rpa

common.sets <- intersect(rownames(emat.old), rownames(emat.new))
common.samples <- intersect(colnames(emat.old), colnames(emat.new))

nd <- emat.new[common.sets, common.samples]
od <- emat.old[common.sets, common.samples]

rs <- sample(prod(dim(nd)), 1e4)
plot(nd[rs], od[rs], pch = ".")

