# ==============================================================================
# Example data for debarcoding
# ------------------------------------------------------------------------------
#' @rdname sample_ff
#' @aliases sample_ff
#' @title 20 sample experiment
#'
#' @description 
#' A flowFrame obtained from a triple-staining experiment.
#' 
#' @details 
#' Following a 6-choose-3 barcoding scheme, mass channels 
#' 102, 104, 105, 106, 108, and 110 were used for labeling
#' such that each of the 20 individual barcodes are positive 
#' for exactly 3 of the 6 possible barcode channels.
#'
#' @docType data
#' @keywords datasets
#' @usage data(sample_ff)
#' @format A flowFrame with 10000 cells and 64 observables.
"sample_ff"

# ==============================================================================
# 6-choose-3 binary barcoding scheme
# ------------------------------------------------------------------------------
#' @rdname sample_key
#' @aliases sample_key
#' @title 6-choose-3 barcoding scheme
#'
#' @description 
#' A data.frame with dimensions 20 x 6.
#' 
#' @details 
#' A 6-choose-3 barcoding scheme, with mass channels 102, 104, 105, 106, 108, 
#' and 110 used for labeling. Each of the 20 individual barcodes are positive 
#' for exactly 3 of the 6 possible barcode channels. For each sample, 
#' the barcoding scheme contains a binary code of length 6, e.g. 111000, 
#' as a unique identifier.
#'
#' @docType data
#' @keywords datasets
#' @usage data(sample_key)
#' @format A data.frame with 20 rows (samples) and 6 columns (barcodes).
"sample_key"

# ==============================================================================
# Example data for compensation 
# ------------------------------------------------------------------------------
#' @rdname ss_beads
#' @aliases ss_beads
#' @title Single-stained beads
#'
#' @description 
#' A flowFrame obtained from a 36ab-panel single-staining experiment.
#' 
#' @details 36ab-panel: Beads were stained for antibodies captured by 
#' mass channels 139, 141 through 156, and 158 through 176, respectively, 
#' and pooled together. Here we sampled a fraction of all recorded events 
#' to demonstrate the debarcoding and compensation work-flow, 
#' at the cost of not necessarily arriving at satisfying results. 
#'
#' @docType data
#' @keywords datasets
#' @usage data(ss_beads)
#' @format A flowFrame with 10000 cells and 61 observables.
"ss_beads"