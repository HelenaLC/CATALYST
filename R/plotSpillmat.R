# ==============================================================================
# Plot spillover matrix heat map
# ------------------------------------------------------------------------------

#' @rdname plotSpillmat
#' @title Spillover matrix heat map
#' 
#' @description 
#' Generates a heat map of the spillover matrix annotated with 
#' estimated spill percentages.
#'
#' @param bc_ms 
#' a vector of numeric masses corresponding to barcode channels.
#' @param SM 
#' spillover matrix returned from \code{computeSpillmat}.
#' @param annotate
#' logical. If TRUE (default), spill percentages are shown inside bins and
#' rows/columns are annotated with the total amount of spill caused/received. 
#' @param palette
#' an optional vector of colors to interpolate.
#' @param out_path 
#' a character string. If specified, outputs will be generated 
#' in this location. Defaults to NULL.
#' @param name_ext 
#' a character string. If specified, will be appended to the plot's name. 
#' Defaults to NULL.
#' 
#' @return plots estimated spill percentages as a heat map. 
#' Colours are ramped to the highest spillover value present
#' 
#' @examples
#' # get single-stained control samples
#' data(ss_exp)
#' 
#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#' 
#' re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' spillMat <- computeSpillmat(x = re)
#' plotSpillmat(bc_ms = bc_ms, SM = spillMat)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom grDevices colorRampPalette
#' @export

# ==============================================================================

plotSpillmat <- function(bc_ms, SM, annotate=TRUE, 
    palette=NULL, out_path=NULL, name_ext=NULL) {
    
    nms <- colnames(SM)
    ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
    bc_cols <- which(ms %in% bc_ms)
    bc_range <- min(bc_cols) : max(bc_cols)
    SM <- make_symetric(SM)[bc_range, bc_range]
    SM <- round(SM*100, 3)
    max <- ceiling(max(SM[SM != 100])/.25)*.25
    suppressWarnings(p <- heatmaply(SM, 
        dendrogram="none", margins=c(100,100,0,0),
        colors=c("white", brewer.pal(9, "Reds")[-1]), 
        limits=c(0, max), na.value="lightgrey", grid_color="white",
        xlab="Receiving", ylab="Emitting",
        colRow=lab_cols,
        label_names=c("Emitting", "Receiving", "Spillover")))

    if (!is.null(out_path)) {
        htmlwidgets::saveWidget(p, file.path(out_path, 
            paste0("SpillMat", name_ext, ".html")))
    } else {
        p
    }
}
