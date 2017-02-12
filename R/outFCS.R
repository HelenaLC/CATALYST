# ==============================================================================
# Plot events
# ------------------------------------------------------------------------------

#' @rdname outFCS
#' @title Write population-wise FCS
#' 
#' @description Writes an FCS file for each barcode population.
#'
#' @param x 
#' a \code{\link{dbFrame}}.
#' @param out_path
#' character string. Specifies in which location 
#' output files are to be generated.
#' @param verbose
# if TRUE (default), a warning is given about populations for which 
#' no FCS files have been generated.
#' 
#' @details 
#' FCS files are not generation for populations with 
#' no or less than 10 event assignments.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' outFCS(x = re, out_path = file.path(tempdir()))
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom flowCore flowFrame write.FCS
#' @export

# ------------------------------------------------------------------------------

setMethod(f="outFCS",       
    signature=signature(x="dbFrame", out_path="character"), 
    definition=function(x, out_path, verbose=TRUE) {
        nms <- rownames(x@bc_key)
        ids <- sort(unique(x@bc_ids))
        skip <- c()
        for (i in ids[ids != 0]) {
            if (sum(x@bc_ids == i) < 10) {
                skip <- c(skip, i) 
                next
            }
            ff <- new("flowFrame", exprs=x@exprs[x@bc_ids == i, ])
            suppressWarnings(
                write.FCS(ff, file.path(out_path, paste0(i, ".fcs"))))
        }
        if (verbose) {
            if (length(skip) > 1) {
                cat("o Samples", paste(skip, collapse=", "), 
                    "contain less than 10 events\n",
                    " (no FCS files have been generated for these barcodes)\n")
            } else if (length(skip) == 1) {
                cat("o Sample", skip, "contains less than 10 events\n",
                    " (no FCS file has been generated for this barcode)\n")
            }
            tmp <- nms[!nms %in% ids]
            if (length(tmp) > 1) {
                cat("o No events assigned to samples", 
                    paste(tmp, collapse=", "),"")
            } else if (length(tmp) == 1) {
                cat("o No events assigned to sample", 
                    paste(tmp, collapse=", "),"")
            }
        }
    })