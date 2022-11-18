#' @rdname plotSpillmat
#' @title Spillover matrix heatmap
#'
#' @description Generates a heatmap of the spillover matrix
#' annotated with estimated spill percentages.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param sm spillover matrix to visualize. If NULL, \code{plotSpillmat} 
#'   will try and access \code{metadata(x)$spillover_matrix}.
#' @param anno logical. If TRUE (default), spill percentages are shown inside 
#'   bins and rows are annotated with the total amount of spill received.
#' @param isotope_list named list. Used to validate the input spillover matrix.
#'   Names should be metals; list elements numeric vectors of their isotopes.
#'   See \code{\link{isotope_list}} for the list of isotopes used by default.
#' @param hm_pal character vector of colors to interpolate.
#' @param anno_col character string specifying 
#'   the color to use for bin annotations.
#'
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#'
#' @return a \code{ggplot2}-object showing estimated spill percentages 
#'   as a heatmap with colors ramped to the highest spillover value present.
#'
#' @examples
#' # get single-stained control samples & construct SCE
#' data(ss_exp)
#' sce <- prepData(ss_exp)
#' 
#' # debarcode single-positive populations
#' bc_ms <- c(139, 141:156, 158:176)
#' sce <- assignPrelim(sce, bc_ms, verbose = FALSE)
#' sce <- applyCutoffs(estCutoffs(sce))
#' 
#' # estimate & visualize spillover matrix
#' sce <- computeSpillmat(sce)
#' plotSpillmat(sce)
#'
#' @import ggplot2 
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @export

plotSpillmat <- function(x, sm = NULL, anno = TRUE, 
    isotope_list = CATALYST::isotope_list,
    hm_pal = c("white", "lightcoral", "red2", "darkred"), anno_col = "black") {
    
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"),
        is.logical(anno), length(anno) == 1,
        !is.null(sm) || !is.null(sm <- metadata(x)$spillover_matrix),
        .check_pal(hm_pal), .check_pal(anno_col, n = 1))

    # check validity of input spillover matrix
    sm <- .check_sm(sm, isotope_list)
    # assure there are no gaps in y-axis channels
    rng <- match(rownames(sm), channels(x))
    chs <- channels(x)[seq(min(rng), max(rng))]
    sm <- .make_symetric(sm)[chs, chs]
    # drop diagonal entries & empty columns
    chs <- rownames(sm)
    sm[cbind(chs, chs)] <- NA
    sm <- sm[, colSums(sm, na.rm = TRUE) != 0]
    
    df <- melt(sm * 100)
    names(df) <- c("from", "to", "spill_val")
    if (hm_pal[1] != "white") hm_pal <- c("white", hm_pal)
    max <- ceiling(max(100 * sm[row(sm) != col(sm)]) / 0.25) * 0.25
    
    p <- ggplot(df, aes_string("to", "from", fill = "spill_val")) +
        geom_tile(col = "lightgrey") +
        scale_fill_gradientn(
            guide = "none", limits = c(0, max),
            colors = hm_pal, na.value = "lightgrey") +
        scale_x_discrete(limits = colnames(sm), expand = c(0, 0)) +
        scale_y_discrete(limits = rev(rownames(sm)), expand = c(0, 0)) +
        coord_fixed() + theme_bw() + theme(
            axis.title = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.border = element_blank(), panel.grid = element_blank())
    if (anno) {
        anno <- sprintf("%.1f", df$spill_val)
        anno[df$spill_val == 0 | is.na(df$spill_val)] <- ""
        p <- p + geom_text(
            aes_string(label = "anno"), 
            col = anno_col, size = 2)
    }
    return(p)
}
