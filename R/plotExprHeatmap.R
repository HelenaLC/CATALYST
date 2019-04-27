#' @rdname plotExprHeatmap
#' @title Plot expression heatmap
#' 
#' @description 
#' Plots median marker expressions across samples 
#' computed on arcsinh-transformed intensities.
#'
#' @param x
#'   a \code{\link{daFrame}}.
#' @param bin_anno 
#'   logical. Specifies whether to display values insinde each bin.
#' @param row_anno 
#'   logical. Should row annotations for each factor 
#'   in \code{metadata(x)} be included?
#' @param palette 
#'   character vector of colors to interpolate. 
#' @param scale 
#'   logical. Specifies whether scaled values should be displayed.
#'   (see below for details)
#' @param draw_freqs 
#'   logical. Specifyies whether to display cell counts and proportions.
#' @param clustering_distance 
#'   a character string that specifies the metric to use in 
#'   \code{\link[stats]{dist}()} for clustering.
#' @param clustering_linkage 
#'   a character string that specifies the linkage to use in
#'   \code{\link[stats]{hclust}()} for clustering.
#' 
#' @details Scaled values corresponds to cofactor arcsinh-transformed 
#' expression values scaled between 0 and 1 using 1% and 99% percentiles as 
#' boundaries. Hierarchical clustering is performed on the unscaled data.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' daf <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotExprHeatmap(daf[, 1:5], draw_freqs=TRUE)
#' 
#' @import ComplexHeatmap SummarizedExperiment
#' @importFrom dplyr funs group_by_ summarize_all select select_if
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata
#' @importFrom scales hue_pal
#' @importFrom stats dist hclust
# ------------------------------------------------------------------------------

setMethod(f="plotExprHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, bin_anno=TRUE, row_anno=TRUE,
        palette=brewer.pal(n=8, name="YlGnBu"), scale=TRUE, draw_freqs=FALSE,  
        clustering_distance="euclidean", clustering_linkage="average") {
        
        # compute medians across samples
        med_exprs <- data.frame(exprs(x), sample_id=sample_ids(x)) %>%
            group_by_(~sample_id) %>% summarize_all(funs(median))
        med_exprs <- data.frame(med_exprs, row.names=1)
        if (scale) 
            med_exprs <- .scale_exprs(med_exprs)
        
        d <- stats::dist(med_exprs, method=clustering_distance)
        row_clustering <- stats::hclust(d, method=clustering_linkage)

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
            Heatmap(matrix=med_exprs, col=hm_cols, name="expression", 
                cell_fun=cell_fun, cluster_rows=row_clustering,
                heatmap_legend_param=list(color_bar="continuous"),
                column_names_gp=gpar(fontsize=8), 
                row_names_gp=gpar(fontsize=8)) 
        }
        if (bin_anno) {
            hm <- hm(cell_fun=function(j, i, x, y, ...)
                grid.text(gp=gpar(fontsize=8),
                    sprintf("%.2f", med_exprs[i, j]), x, y))
        } else {
            hm <- hm(cell_fun=function(...) NULL)
        }
        
        if (row_anno) {
            md <- metadata(x)$experiment_info
            m <- match(rownames(med_exprs), md$sample_id)
            df <- select(md[m, ], -"sample_id")
            df <- select_if(df, is.factor)
            row_anno <- .anno_factors(df, "row")
            row_anno + hm + freq_bars + freq_anno
        } else {
            hm + freq_bars + freq_anno
        }
    }
)
