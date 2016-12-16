# load single staining experiment 
path <- "/Users/HLC/Dropbox/spillover/corrected bead based compensation/Exp3"
fcs <- file.path(path, "160805_Exp3_beads-before.fcs")
library(flowCore)
ss_exp  <- read.FCS(fcs)

# specify mass channels stained for
bc_ms <- c(139, 141:156, 158:176)

# ------------------------------------------------------------------------------
# NOTE: may want to specify option out_path where indicated
# ------------------------------------------------------------------------------

# assign preliminary IDs
library(CATALYST)
re <- assignPrelim(x = ss_exp, y = bc_ms) 
plotEvents(x = re, which_bc = "all", n_events = 100, name_ext = "prelim") # out_path = ...

# estimate population separation cutoffs
re <- estCutoffs(x = re)  
plotYields(x = re, which_bc = c(0, bc_ms)) # out_path = ...

# optionally, choose cutoffs manually
#sep_cutoffs(re) <- ...

# apply deconvolution paramters
re <- applyCutoffs(x = re, mhl_cutoff = 10)
plotEvents(x = re, which_bc = "all", n_events = 100, name_ext = "_final") # out_path = ...

# write popultion-wise FCS
#outFCS(x = re, out_path = ...)

# estimate optimal trim value
trimVal <- estTrim(x = re, .06, .14, .02) # out_path = ...

# estimate compensation matrix
spillMat <- computeSpillmat(x = re, method = "mean", trim = .08) 

# plot spillover matrix heat map
plotSpillmat(bc_ms = bc_ms, SM = spillMat) # out_path = ...

# compensate single staining 
comped_ss_beads <- ss_beads
exprs(comped_ss_beads) <- exprs(ss_beads) %*% compMat
write.FCS(comped_ss_beads, filename = paste0(gsub(".fcs", "", fcs), "_comped.fcs"))

# compensate all FCS files in a folder
# !!!need to specify path to folder!!!
compCytof(path, spillMat)
    
