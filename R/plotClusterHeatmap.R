#' @rdname plotClusterHeatmap
#' @title Plot cluster heatmap
#' 
#' @description Plots expression & relative cluster abundances 
#' heatmaps summarizing a clustering and/or metaclustering.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param hm2 character string. Specifies the right-hand side heatmap. 
#'   One of: \itemize{
#'   \item{\code{"abundances"}: cluster frequencies across samples}
#'   \item{\code{"state"}: median cell state marker expressions 
#'     across clusters (analogous to the left-hand side heatmap)}
#'   \item{a character string/vector corresponding to one/multiple marker(s): 
#'     aggregated marker expressions across samples and clusters}}
#' @param k character string specifying the clustering 
#'   across which median marker expressions should be computed.
#' @param m character string specifying the metaclustering to be shown. 
#'   (This is for display only and will not effect any computations!) 
#' @param fun character string specifying 
#'   the function to use as summary statistic.
#' @param cluster_anno logical specifying if clusters should be annotated.
#' @param split_by deprecated.
#' @param scale logical specifying whether scaled values should be plotted.
#' @param draw_dend logical specifying if the row dendrogram should be drawn.
#' @param draw_freqs logical specifying 
#'   whether to display cell counts and proportions.
#' @param palette character vector of colors to interpolate. 
#'   
#' @details Scaled values corresponds to cofactor arcsinh-transformed 
#' expression values scaled between 0 and 1 using 1% and 99% percentiles as 
#' boundaries. Hierarchical clustering is performed on the unscaled data.
#' 
#' In its 1st panel, \code{plotClusterHeatmap} will display
#' median (scaled, arcsinh-transformed) cell-type marker expressions (across all samples).
#' Depending on argument \code{hm2}, the 2nd panel will contain one of:
#' \itemize{
#' \item{relataive cluster abundances by sample}
#' \item{median (scaled, arcsinh-transformed) 
#'   cell-state marker expressions (across all samples)}
#' \item{median (scaled, arcsinh-transformed) 
#'   cell-state marker expressions by sample}
#' }
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @export

plotClusterHeatmap <- function(x, hm2 = NULL, 
    k = "meta20", m = NULL, fun = c("median", "mean"), 
    cluster_anno = TRUE, split_by = NULL, scale = TRUE, 
    draw_dend = TRUE, draw_freqs = FALSE, 
    palette = rev(brewer.pal(11, "RdYlBu"))) {
    
    .Deprecated(
        new = "plotMultiHeatmap",
        old = "plotClusterHeatmap",
        msg = paste(sep = "\n",
            "'plotClusterHeatmap' is deprecated; instead, please use",
            " o 'plotExprHeatmap' for aggregated expression heatmaps",
            " o 'plotFreqHeatmap' for cluster frequency heatmaps",
            " o 'plotMultiHeatmap' to combine multiple heatmaps"))
    
    if (is.null(hm2)) {
        plotExprHeatmap(x, features = "type", 
            by = "cluster_id", k = k, m = m,
            assay = "exprs", fun = match.arg(fun),
            scale = "first", q = 0.01,
            row_anno = cluster_anno, col_anno = FALSE,
            row_clust = TRUE, col_clust = FALSE,
            row_dend = TRUE, col_dend = FALSE,
            bars = draw_freqs, perc = draw_freqs, 
            hm_pal = palette)
    } else {
        plotMultiHeatmap(x, 
            hm1 = "type", hm2 = hm2, 
            k = k, m = m, 
            assay = "exprs", fun = match.arg(fun), 
            scale = ifelse(scale, "first", "never"), 
            q = 0.01, normalize = FALSE,
            row_anno = cluster_anno, col_anno = FALSE, 
            row_clust = TRUE, col_clust = FALSE, 
            row_dend = TRUE, col_dend = FALSE, 
            bars = draw_freqs, perc = draw_freqs, 
            hm1_pal = palette)
    }
}
