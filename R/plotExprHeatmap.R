#' @rdname plotExprHeatmap
#' @title Plot expression heatmap
#' 
#' @description 
#' Heatmap of median marker expressions by sample, include annotation 
#' of cell metadata factors as well as relative and absolute cell counts.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param scale logical specifying whether expressions should be scaled
#'   using lower (1\%) and upper (99\%) quantiles as boundaries;
#'   hierarchical clustering is performed on unscaled data, 
#'   regardless of whether \code{scale = TRUE} or \code{FALSE}.
#' @param row_anno logical. Should row annotations for 
#'   each non-numeric cell metadata factor be included?
#' @param row_clust,col_clust logical specifying 
#'   whether rows/columns (samples/features) should be 
#'   hierarchically clustered and re-ordered accordingly.
#' @param row_dend,col_dend logical specifying whether to include the
#'   row/column dendrogram for the hierarchical clustering of samples/markers.
#' @param bin_anno logical. Specifies whether to display values inside bins.
#' @param draw_freqs logical specifying whether to display
#'   a barplot of cell counts labeled with proportions 
#'   for each sample as a right-hand side row annotation.
#' @param hm_pal character vector of colors to interpolate for the heatmap. 
#' @param distance character string specifying the distance metric 
#'   to use in \code{\link[stats]{dist}} for hierarchical clustering. 
#' @param linkage character string specifying the agglomeration method 
#'   to use in \code{\link[stats]{hclust}} for hierarchical clustering. 
#' 
#' @return a \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Crowell HL, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # turn everything on
#' plotExprHeatmap(sce, 
#'   bin_anno = TRUE,
#'   draw_freqs = TRUE)
#'   
#' # turn everything off
#' plotExprHeatmap(sce,
#'   row_anno = FALSE,
#'   row_dend = FALSE,
#'   col_dend = FALSE)
#' 
#' @import ComplexHeatmap SummarizedExperiment
#' @importFrom dplyr select select_if
#' @importFrom grid grid.text
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay
#' @export

plotExprHeatmap <- function(x, scale = TRUE, row_anno = TRUE, 
    row_clust = TRUE, col_clust = TRUE, row_dend = TRUE, col_dend = TRUE, 
    draw_freqs = FALSE,  bin_anno = FALSE, 
    hm_pal = brewer.pal(9, "YlGnBu"), 
    distance = c(
        "euclidean", "maximum", "manhattan", 
        "canberra", "binary", "minkowski"), 
    linkage = c(
        "average", "ward.D", "single", "complete", 
        "mcquitty", "median", "centroid", "ward.D2")) {
    
    # check validity of input arguments
    .check_sce(x)
    .check_colors(hm_pal)
    stopifnot(
        is.logical(scale), length(scale) == 1,
        is.logical(bin_anno), length(bin_anno) == 1,
        is.logical(row_anno), length(row_anno) == 1,
        is.logical(row_dend), length(row_dend) == 1,
        is.logical(col_dend), length(col_dend) == 1,
        is.logical(row_clust), length(row_clust) == 1,
        is.logical(col_clust), length(col_clust) == 1,
        is.logical(draw_freqs), length(draw_freqs) == 1)
    distance <- match.arg(distance)
    linkage <- match.arg(linkage)
    
    # compute medians across samples
    z <- t(.agg(x, "sample_id"))
    # (optionally) cluster rows (markers)
    if (row_clust) {
        d <- dist(z, method = distance)
        row_clust <- hclust(d, method = linkage) 
    }
    # (optionally) do 0-1 scaling for each marker
    if (scale) 
        z <- .scale_exprs(z, 2)
    
    # left-hand side heatmap annotation:
    # non-numeric cell metadata variables
    cd <- data.frame(colData(x))
    cd <- select_if(cd, is.factor)
    cols_keep <- setdiff(colnames(cd), 
        c("cluster_id", "sample_id"))
    if (row_anno && length(cols_keep) > 0) {
        m <- match(rownames(z), x$sample_id)
        cd <- cd[m, cols_keep, drop = FALSE]
        left_anno <- .anno_factors(cd, "row")
    } else left_anno <- NULL
    
    # right-hand side heatmap annotation:
    # labeled barplot of event counts by sample
    freq_bars <- freq_anno <- NULL
    if (draw_freqs) {
        ns <- tabulate(x$sample_id)
        fq <- round(ns/sum(ns)*100, 2)
        txt <- sprintf("%s (%s%%)", rownames(z), fq)
        right_anno <- rowAnnotation(
            "n_cells" = row_anno_barplot(
                x = ns, width = max_text_width(txt),
                gp = gpar(fill = "grey", col = "white"),
                border = FALSE, axis = TRUE, bar_width = 0.8),
            "foo" = row_anno_text(x = txt, location = 0.5, just = "centre"))
    } else right_anno <- NULL
    
    if (bin_anno) {
        cell_fun <- function(j, i, x, y, ...) 
            grid.text(
                gp = gpar(fontsize = 8), 
                sprintf("%.2f", z[i, j]), x, y)
    } else cell_fun <- NULL
    
    Heatmap(
        matrix = z, 
        name = paste0("scaled\n"[scale], "expression"),
        col = colorRampPalette(hm_pal)(100), 
        cell_fun = cell_fun, 
        cluster_rows = row_clust,
        cluster_columns = col_clust,
        show_row_dend = row_dend,
        show_column_dend = col_dend,
        show_row_names = is.null(right_anno),
        row_names_side = ifelse(row_anno || row_dend, "right", "left"),
        left_annotation = left_anno,
        right_annotation = right_anno)
}
