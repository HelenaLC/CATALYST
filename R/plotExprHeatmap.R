#' @rdname plotExprHeatmap
#' @title Plot expression heatmap
#' 
#' @description 
#' Plots median marker expressions across samples 
#' computed on arcsinh-transformed intensities.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param bin_anno logical. Specifies whether to display values inside bins.
#' @param row_anno logical. Should row annotations for each factor 
#'   in \code{metadata(x)$experiment_info} be included?
#' @param palette character vector of colors to interpolate. 
#' @param scale logical. Should scaled values be displayed? (see details)
#' @param draw_freqs logical. Should cell counts and proportions be displayed?
#' @param clustering_distance character string that specifies 
#'   the metric to use in \code{\link[stats]{dist}} for clustering.
#' @param clustering_linkage character string that specifies 
#'   the linkage to use in \code{\link[stats]{hclust}} for clustering.
#' 
#' @details Scaled values corresponds to cofactor arcsinh-transformed 
#' expression values scaled between 0 and 1 using 1% and 99% percentiles as 
#' boundaries. Hierarchical clustering is performed on the unscaled data.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
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
#' plotExprHeatmap(sce, draw_freqs=TRUE)
#' 
#' @import ComplexHeatmap SummarizedExperiment
#' @importFrom dplyr select select_if
#' @importFrom grid grid.text
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata
#' @importFrom scales hue_pal
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay
#' @export

plotExprHeatmap <- function(x, bin_anno=TRUE, row_anno=TRUE,
    palette=brewer.pal(n=8, name="YlGnBu"), scale=TRUE, draw_freqs=FALSE,  
    clustering_distance="euclidean", clustering_linkage="average") {
    
    .check_sce(x)
    
    # compute medians across samples
    ms <- t(.agg(x, "sample_id"))
    d <- dist(ms, method=clustering_distance)
    row_clustering <- hclust(d, method=clustering_linkage)
    if (scale) ms <- .scale_exprs(ms, 2)
    
    # barplots of event counts
    freq_bars <- freq_anno <- NULL
    if (draw_freqs) {
        counts <- as.numeric(n_cells(x))
        freqs <- round(counts/sum(counts)*100, 2)
        freq_bars <- rowAnnotation(width=unit(2, "cm"), 
            "n_cells"=row_anno_barplot(x=counts, border=FALSE, axis=TRUE,
                gp=gpar(fill="grey50", col="white"), bar_with=.8))
        labs <- paste0(counts, " (", freqs, "%)")
        freq_anno <- rowAnnotation(
            text=row_anno_text(labs), 
            width=max_text_width(labs))
    } 
    
    # heatmap of medians across antigens and samples
    hm_cols <- colorRampPalette(palette)(100)
    hm <- function (cell_fun) { 
        Heatmap(matrix=ms, col=hm_cols, name="expression", 
            cell_fun=cell_fun, cluster_rows=row_clustering,
            heatmap_legend_param=list(color_bar="continuous"),
            column_names_gp=gpar(fontsize=8), 
            row_names_gp=gpar(fontsize=8)) 
    }
    if (bin_anno) {
        hm <- hm(cell_fun=function(j, i, x, y, ...)
            grid.text(gp=gpar(fontsize=8),
                sprintf("%.2f", ms[i, j]), x, y))
    } else {
        hm <- hm(cell_fun=function(...) NULL)
    }
    
    if (row_anno) {
        md <- metadata(x)$experiment_info
        m <- match(rownames(ms), md$sample_id)
        df <- select(md[m, ], -"sample_id")
        df <- select_if(df, is.factor)
        row_anno <- .anno_factors(df, "row")
        row_anno + hm + freq_bars + freq_anno
    } else {
        hm + freq_bars + freq_anno
    }
}