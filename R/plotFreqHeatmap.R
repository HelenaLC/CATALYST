# ==============================================================================
# Heatmap for differental abundance & state analysis
# ------------------------------------------------------------------------------
#' @rdname plotFreqHeatmap
#' @title Cluster frequency heatmap
#' @description Heatmap of relative cluster abundances (frequencies) by sample.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k character string specifying the clustering to use;
#'   valid values are \code{names(cluster_codes(x))}.
#'   Cell counts will be computed across these cluster IDs.
#' @param m character string specifying a metaclustering 
#'   to include as an annotation when \code{row_anno = TRUE}.
#' @param normalize logical specifying whether to Z-score normalize. 
#' @param row_anno,col_anno logical specifying whether to 
#'   include row/column annotations for clusters/samples; 
#'   for \code{col_anno}, this can be a character vector specifying 
#'   a subset of \code{names(colData(x))} to be included.
#' @param row_clust,col_clust logical specifying 
#'   whether rows/columns (clusters/samples) should be 
#'   hierarchically clustered and re-ordered accordingly.
#' @param row_dend,col_dend logical specifying 
#'   whether to include row/column dendrograms.
#' @param bars logical specifying whether to include a barplot 
#'   of cell counts per cluster as a right-hand side row annotation.
#' @param perc logical specifying whether to display 
#'   percentage labels next to bars when \code{bars = TRUE}.
#' @param hm_pal character vector of colors to interpolate for the heatmap. 
#' @param k_pal,m_pal character vector of colors
#'   to use for cluster and merging row annotations.
#'   If less than \code{nlevels(cluster_ids(x, k/m))} 
#'   values are supplied, colors will be interpolated 
#'   via \code{\link[grDevices:colorRamp]{colorRampPalette}}.
#' 
#' @return a \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @seealso 
#' \code{\link{plotAbundances}}, 
#' \code{\link{plotExprHeatmap}}, 
#' \code{\link{plotMultiHeatmap}},
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' # complete 
#' plotFreqHeatmap(sce, k = "meta12", m = "meta8")
#' 
#' # minimal 
#' plotFreqHeatmap(sce, k = "meta10",
#'   normalize = FALSE, bars = FALSE,
#'   row_anno = FALSE, col_anno = FALSE,
#'   row_clust = FALSE, col_clust = FALSE)
#'   
#' # customize colors & annotations
#' plotFreqHeatmap(sce, 
#'   k = "meta7", m = "meta4",
#'   col_anno = "condition",
#'   hm_pal = c("navy", "grey95", "gold"),
#'   k_pal = hcl.colors(7, "Set 2"),
#'   m_pal = hcl.colors(4, "Dark 3"))
#' 
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom Matrix colSums
#' @importFrom RColorBrewer brewer.pal
#' @export

plotFreqHeatmap <- function(x, 
    k = "meta20", m = NULL, normalize = TRUE, 
    row_anno = TRUE, col_anno = TRUE,
    row_clust = TRUE, col_clust = TRUE, 
    row_dend = TRUE, col_dend = TRUE,
    bars = TRUE, perc = FALSE,
    hm_pal = rev(brewer.pal(11, "RdBu")),
    k_pal = .cluster_cols, m_pal = k_pal) {
    
    # check validity of input arguments
    args <- as.list(environment())
    .check_args_plotFreqHeatmap(args)

    # compute cell counts & frequencies by cluster-sample
    x$cluster_id <- cluster_ids(x, k)
    ns <- table(x$cluster_id, x$sample_id)
    fq <- prop.table(ns, 2)
    y <- as.matrix(unclass(fq))
    
    # do Z-score normalization by cluster & across samples
    if (normalize) y <- .z_normalize(asin(sqrt(y)))
    
    # left-hand side row annotation of (meta)clustering
    if (!isFALSE(row_anno)) {
        left_anno <- .anno_clusters(x, k, m, k_pal, m_pal)
    } else left_anno <- NULL
    # column annotation of cell metadata variables
    if (!isFALSE(col_anno)) {
        sids <- levels(droplevels(factor(x$sample_id)))
        top_anno <- .anno_factors(x, sids, col_anno, "colum") 
    } else top_anno <- NULL
 
    # right-hand side row annotation of cell counts & frequencies by cluster
    if (bars) {
        right_anno <- .anno_counts(x$cluster_id, perc)
    } else right_anno <- NULL
    
    Heatmap(
        matrix = y, 
        name = paste0("normalized\n"[normalize], "frequency"),
        col = hm_pal,
        na_col = "lightgrey", 
        rect_gp = gpar(col = "white"),  
        column_title = "sample_id",
        column_title_side = "bottom",
        cluster_rows = row_clust,
        cluster_columns = col_clust,
        show_row_dend = row_dend,
        show_column_dend = col_dend,
        show_row_names = is.null(left_anno),
        row_names_side = "left",
        top_annotation = top_anno,
        left_annotation = left_anno,
        right_annotation = right_anno)
}
