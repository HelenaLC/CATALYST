library(CATALYST)
library(flowCore)

path <- "/Users/HLC/Dropbox/spillover/corrected bead based compensation/Exp3"
beads_ff  <- read.FCS(file.path(path, "160805_Exp3_beads-before.fcs"))

bc_ms <- c(139, 141:156, 158:176)
re <- assignPrelim(x = beads_ff, y = bc_ms) 
re <- estCutoffs(x = re)   
re <- applyCutoffs(x = re, mhl_cutoff = 10)

trimVal <- estTrim(x = re, .06, .14, .02)
compMat <- computeCompmat(x = re, method = "mean", trim = trimVal) 
plotSpillmat(bc_ms = bc_ms, CM = compMat)
plotScatter(x = re, CM = compMat)