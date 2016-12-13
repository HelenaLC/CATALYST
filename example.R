library(CATALYST)
library(flowCore)

# ==============================================================================
# compute debarcoding in one line
# ------------------------------------------------------------------------------
data(ss_beads)
debarcode(x = ss_beads, y = bc_ms, out_path = "~/Desktop")
# ------------------------------------------------------------------------------

path <- "/Users/HLC/Dropbox/spillover/corrected bead based compensation/Exp3"
beads_ff  <- read.FCS(file.path(path, "160805_Exp3_beads-before.fcs"))

# specify barcode key:
# - binary table w/ masses as column names OR
# - vector of masses for single-staining experiment
bc_ms <- c(139, 141:156, 158:176)

# assign preliminary IDs
re <- assignPrelim(x = beads_ff, y = bc_ms) 
plotEvents(x = re, which_bc = 172, n_events = 20)

# estimate population separation cutoffs
re <- estCutoffs(x = re)  
plotYields(x = re, which_bc = c(0, 139))

# apply deconvolution paramters
re <- applyCutoffs(x = re, mhl_cutoff = 10)
plotEvents(x = re, which_bc = 172, n_events = 20)

# write popultion-wise FCS
#outFCS(x = re, out_path = ...)

# estimate optimal trim value
trimVal <- estTrim(x = re, .06, .14, .02)

# estimate compensation matrix
compMat <- computeCompmat(x = re, method = "mean", trim = trimVal) 

# before-vs-after compensation at a glance
plotSpillmat(bc_ms = bc_ms, CM = compMat)
