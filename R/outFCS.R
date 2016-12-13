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
#' 
#' @examples
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom flowCore flowFrame write.FCS
#' @export

# ------------------------------------------------------------------------------

setMethod(f="outFCS",       
    signature=signature(x="dbFrame", out_path="character"), 
    definition=function(x, out_path) {
        nms <- rownames(x@bc_key)
        ids <- sort(unique(x@bc_ids))
        for (i in ids[ids != 0]) {
            if (sum(x@bc_ids == i) < 10) {
                cat("Sample", i, "contains less than 10 events;",
                    "no FCS file has been generated.\n")
                next
            }
            ff <- new("flowFrame", exprs=x@exprs[x@bc_ids == i, ])
            write.FCS(ff, file.path(out_path, paste0(i, ".fcs")))
        }
        cat("No events assigned to sample(s)", 
            paste(nms[!nms %in% ids], collapse=", "))
    })
        