# ================================================================================
# Plot events
# --------------------------------------------------------------------------------

#' @rdname outFCS
#' @title Write population-wise FCS
#' 
#' @description Writes an FCS file for each barcode population.
#'
#' @param x        a \code{\link{dbFrame}}.
#' @param out_path character string. Specifies in which location output files are to be generated.
#' @param names    optional vector of strings specifying output file names. 
#'                 By default, each FCS will be given the assigned ID as name.
#' 
#' @examples
#' data(ss_beads)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_beads, y = bc_ms)
#' #outFCS(x = re, out_path = ...)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom flowCore flowFrame write.FCS
#' @export

# --------------------------------------------------------------------------------

setMethod(f="outFCS",       
    signature=signature(x="dbFrame"), 
    definition=function(x, out_path, names=NULL) {
        if (is.null(names)) 
            names <- rownames(x@bc_key)
        ids <- sort(unique(x@bc_ids))
        ids <- ids[ids != 0]
        for (i in ids) {
            ff <- new("flowFrame", exprs=x@exprs[x@bc_ids == i, ])
            write.FCS(ff, file.path(out_path, paste0(names[ids == i], ".fcs")))
        }
    })