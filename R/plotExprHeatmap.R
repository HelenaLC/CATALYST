#' @rdname plotExprHeatmap
#' @title Plot expression heatmap
#' 
#' @description 
#' Plots median marker expressions across samples 
#' computed on arcsinh-transformed intensities.
#'
#' @param x
#'   a \code{\link{daFrame}}.
#' @param anno 
#'   logical. Specifies whether to display values insinde each bin.
#' @param color_by 
#'   character string. Specifies the row annotation.
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
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotExprHeatmap(re[, 1:5], draw_freqs=TRUE)
#' 
#' @import ComplexHeatmap SummarizedExperiment
#' @importFrom dplyr funs group_by_ summarize_all
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata
#' @importFrom scales hue_pal
#' @importFrom stats dist hclust
# ------------------------------------------------------------------------------

setMethod(f="plotExprHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, anno=TRUE, color_by=NULL,
        palette=brewer.pal(n=8, name="YlGnBu"), scale=TRUE, draw_freqs=FALSE,  
        clustering_distance="euclidean", clustering_linkage="average") {

        # validity check
        valid_opts <- colnames(rowData(x))
        if (!is.null(color_by) && !color_by %in% valid_opts)
            stop("Invalid argument 'color_by'.\nShould be one of ", 
                paste(dQuote(valid_opts), collapse=", "), " or NULL.")
        
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
            counts <- as.numeric(metadata(x)$n_cells)
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
        if (anno) {
            hm <- hm(cell_fun=function(j, i, x, y, ...)
                grid.text(gp=gpar(fontsize=8),
                    sprintf("%.2f", med_exprs[i, j]), x, y))
        } else {
            hm <- hm(cell_fun=function(...) NULL)
        }
        
        if (!is.null(color_by)) {
            md <- metadata(x)$experiment_info
            m <- match(rownames(med_exprs), md$sample_id)
            row_anno <- data.frame(md[m, color_by], row.names=md$sample_id[m])
            names(row_anno) <- color_by
            row_anno <- Heatmap(matrix=row_anno, name=color_by,
                col=scales::hue_pal()(nlevels(md[, color_by])),
                cluster_rows=row_clustering, show_row_names=FALSE, 
                rect_gp=gpar(col="white"), width=unit(.5, "cm"))
            row_anno + hm + freq_bars + freq_anno
        } else {
            hm + freq_bars + freq_anno
        }
    }
)
