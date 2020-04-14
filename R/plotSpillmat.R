#' @rdname plotSpillmat
#' @title Spillover matrix heatmap
#'
#' @description Generates a heatmap of the spillover matrix
#' annotated with estimated spill percentages.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param sm spillover matrix to visualize. If NULL, \code{plotSpillmat} 
#'   will try and access \code{metadata(x)$spillover_matrix}.
#' @param out_path character string. If specified, outputs will be generated here.
#' @param name_ext character string. If specified, will be appended to the plot's name. 
#' @param anno logical. If TRUE (default), spill percentages are shown inside 
#'   bins and rows are annotated with the total amount of spill received.
#' @param plotly logical. Should an interactive plot be rendered?
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
#' @importFrom htmltools save_html
#' @importFrom methods is
#' @importFrom plotly ggplotly layout
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @export

plotSpillmat <- function(x, sm = NULL, out_path = NULL, name_ext = NULL,
    anno = TRUE, plotly = FALSE, isotope_list = CATALYST::isotope_list,
    hm_pal = c("white", "lightcoral", "red2", "darkred"), anno_col = "black") {
    
    # check validity of input arguments
    .check_colors(hm_pal)
    .check_colors(anno_col, n = 1)
    stopifnot(is(x, "SingleCellExperiment"),
        !is.null(sm) || !is.null(sm <- metadata(x)$spillover_matrix),
        is.null(name_ext) || is.character(name_ext) && length(name_ext) == 1,
        is.null(out_path) || dir.exists(out_path),
        is.logical(anno), length(anno) == 1,
        is.logical(plotly), length(plotly) == 1)

    # check validity of input spillover matrix
    sm <- .check_sm(sm, isotope_list)
    ms <- .get_ms_from_chs(chs <- colnames(sm))
    chs <- rowData(x)$channel_name
    bc_chs <- chs[rowData(x)$is_bc]
    bc_idx <- which(bc_chs %in% colnames(sm))
    bc_rng <- seq(min(bc_idx), max(bc_idx))
    sm <- .make_symetric(sm)[bc_rng, bc_rng]
    
    df <- melt(sm * 100)
    colnames(df) <- c("emitting", "receiving", "spill")
    df$spillover <- paste0(sprintf("%2.3f", df$spill), "%")
    max <- ceiling(max(100 * sm[row(sm) != col(sm)]) / 0.25) * 0.25
    total <- paste0(sprintf("%2.2f", colSums(sm) * 100 - 100, "%"))
    
    labs <- chs[bc_rng]
    ex <- !labs %in% chs[bc_idx]
    lab_cols <- rep("black", nrow(sm))
    lab_cols[ex] <- "grey"
    total[ex] <- NA
    
    if (hm_pal[1] != "white")
        hm_pal <- c("white", hm_pal)
    p <- ggplot(df, 
        aes_string("receiving", "emitting", group = "spillover")) +
        geom_tile(aes_string(fill = "spill"), col = "lightgrey") +
        scale_fill_gradientn(
            guide = FALSE, limits = c(0, max),
            colors = hm_pal, na.value = "lightgrey") +
        scale_x_discrete(limits = colnames(sm), expand = c(0, 0)) +
        scale_y_discrete(limits = rev(rownames(sm)), expand = c(0, 0)) +
        coord_fixed() + theme_bw() + theme(
            axis.title = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.border = element_blank(), panel.grid = element_blank())
    if (anno) {
        anno <- sprintf("%.1f", df$spill)
        anno[df$spill == 0 | df$spill == 100] <- ""
        p <- p + geom_text(
            aes_string(label = "anno"), col = anno_col,
            size = ifelse(is.null(out_path), 2, 3))
    }
    if (plotly)
        p <- ggplotly(p, width = 720, height = 720,
            tooltip = c("emitting", "receiving", "spillover"))
    
    if (is.null(out_path)) {
        p 
    } else {
        ext <- ifelse(plotly, ".html", ".pdf")
        fn <- paste0("spillover_matrix", name_ext, ext)
        fn <- file.path(out_path, fn) 
        if (is(p, "plotly")) {
            save_html(p, fn)
        } else {
            ggsave(fn, p, width = 7, height = 7)
        }
    } 
}
