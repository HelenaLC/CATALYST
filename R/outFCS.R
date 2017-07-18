# ==============================================================================
# Plot events
# ------------------------------------------------------------------------------
#' @rdname outFCS
#' @title Write population-wise FCS files
#' @description Writes an FCS file for each sample from a dbFrame.
#'
#' @param x 
#' a \code{\link{dbFrame}}.
#' @param y
#' a \code{\link{flowFrame}} containing the original measurement and meta data.
#' @param out_path
#' character string. Specifies in which location 
#' output files are to be generated.
#' @param out_nms
#' an optional character string. Either the name of a 2 column CSV table 
#' with sample IDs and desired output file names, or a vector of length 
#' \code{nrow(bc_key(x))} ordered as the samples in the barcoding scheme. 
#' If NULL (default), sample IDs will be used as file names.
#' @param verbose 
#' if TRUE (default), a warning is given about populations 
#' for which no FCS files have been generated.
#' 
#' @details 
#' Creates a separate FCS file for each barcode population. If \code{out_nms} 
#' is NULL (the default), files will be named after the barcode population's ID 
#' in the \code{bc_key} slot of the input \code{\link{dbFrame}}; 
#' unassigned events will be written to "unassigned.fcs", and no output 
#' is generated for populations with less than 10 event assignments.
#' 
#' @return a character of the output path.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' outFCS(x = re, y = sample_ff)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom flowCore flowFrame write.FCS
#' @importFrom utils read.csv
#' @export
# ------------------------------------------------------------------------------

setMethod(f="outFCS", 
    signature=signature(x="dbFrame", y="flowFrame"),
    definition=function(x, y, out_path=tempdir(), out_nms=NULL, verbose=TRUE) {
        stopifnot(is.character(out_path), length(out_path) == 1L)
        smpl_nms <- rownames(bc_key(x))
        if (is.null(out_nms)) {
            out_nms <- rownames(bc_key(x))
        } else if (is.character(out_nms)) {
            if (is.null(dim(out_nms))) {
                if (length(out_nms) != nrow(bc_key(x)))
                    stop("Only ", length(out_nms), " file name(s)",
                        " provided but ", nrow(bc_key(x)), " needed.")
                out_nms <- out_nms
            } else {
                nms_tbl <- utils::read.csv(out_nms, header=FALSE)
                if (nrow(nms_tbl) != nrow(bc_key(x)))
                    stop("Only ", nrow(nms_tbl), " file name(s)", 
                        " provided but ", nrow(bc_key(x)), " needed.")
                if (sum(smpl_nms %in% nms_tbl[, 1]) != nrow(bc_key(x))) 
                    stop("Couldn't find a file name for all samples.",
                        "\nPlease make sure all sample IDs occur ", 
                        "in the provided naming scheme.")
                out_nms <- paste0(nms_tbl[,2], "_", smpl_nms)
            }
        }
        ids <- sort(unique(bc_ids(x)))
        inds <- match(bc_ids(x), ids)
        for (i in seq_along(ids)) {
            id <- ids[i]
            if (id == "0") {
                nm <- "Unassigned"
            } else {
                nm <- out_nms[smpl_nms == id]
            }
            suppressWarnings(flowCore::write.FCS(x=y[inds == i, ], 
                filename=file.path(out_path, paste0(nm, ".fcs"))))
        }
        if (verbose) {
            empty <- smpl_nms[!smpl_nms %in% ids]
            if (length(empty) > 1)
                message("o No events assigned to sample(s) ", 
                    paste(empty, collapse=", "))
        }
        message("*** ", length(ids), " FCS files created in")
        out_path
    })