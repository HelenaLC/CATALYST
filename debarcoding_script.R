# load data and barcoding scheme
library(flowCore)
ff <- read.FCS("/Users/HLC/Dropbox/spillover/debarcoder/20161212_SC_titration_1.FCS")
key <- read.csv("/Users/HLC/Dropbox/spillover/debarcoder/bc_key.csv", 
    check.names = FALSE, row.names = 1)

# assign preliminary IDs
library(CATALYST)
re <- assignPrelim(x = ff, y = key) 
plotEvents(x = re, which = "all", n_events = 100)

# estimate population separation cutoffs
re <- estCutoffs(x = re)  
plotYields(x = re, which = c(0, sort(unique(bc_ids(re)))))

# optionally, choose cutoffs manually
#sep_cutoffs(re) <- ...

# apply deconvolution parameters
re <- applyCutoffs(x = re, mhl_cutoff = 25)
plotMahal(x = re, which = "B3")
plotEvents(x = re, which = "all", n_events = 100)

# write popultion-wise FCS
#outFCS(x = re, out_path = ...)
