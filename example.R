library(CATALYST)
library(flowCore)

# ==============================================================================
# compute debarcoding in one line
# ------------------------------------------------------------------------------
ff <- read.FCS("/Users/HLC/Dropbox/spillover/debarcoder/20161212_SC_titration_1.FCS")
key <- read.csv("/Users/HLC/Dropbox/spillover/debarcoder/bc_key.csv", 
    check.names = FALSE, row.names = 1)
debarcode(x = ff, y = key, out_path = "~/Desktop")


# ==============================================================================
# debarcode and compensate step-by-step
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
