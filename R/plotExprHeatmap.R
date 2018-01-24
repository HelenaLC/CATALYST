# ==============================================================================
# Heatmap of median marker expressions across samples
# ------------------------------------------------------------------------------
#' @rdname plotExprHeatmap
#' @title Median marker expressions across samples
#' 
#' @description Plots median marker expressions across samples
#' computed on arcsinh-transformed intensities.
#'
#' @param x expression matrix.
#' @param anno logical. Specifies whether to display values insinde each bin.
#' @param palette a character string specifying the colors to interpolate.
#' @param out_path a character string. If specified, 
#' output will be generated in this location. Defaults to NULL.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotExprHeatmap(re)
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import ComplexHeatmap SummarizedExperiment
#' @importFrom dplyr funs group_by summarize_all
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats dist hclust
#' @export
# ==============================================================================

setMethod(f="plotExprHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, anno=TRUE, 
        palette=brewer.pal(n=8, name="YlGnBu"), out_path=NULL) {
        
        es <- exprs(x)
        md <- metadata(x)[[1]]
        
        # compute medians across samples
        med_exprs <- data.frame(es, sample_id=sample_ids(x)) %>%
            group_by(sample_id) %>% summarize_all(funs(median))
        med_exprs <- data.frame(med_exprs, row.names=1)
        
        d <- stats::dist(med_exprs, method="euclidean")
        row_clustering <- stats::hclust(d, method="average")
        
        # row annotations
        m <- match(rownames(med_exprs), md$sample_id)
        cond_labs <- data.frame(condtion=md$condition[m],
            row.names=rownames(med_exprs))
        levels(cond_labs) <- levels(md$condition)
        row_anno <- Heatmap(
            cond_labs, c("#F8766D", "#00BFC4"), "condition", 
            cluster_rows=row_clustering, show_row_names=FALSE, 
            rect_gp=gpar(col="white"), width=unit(.5, "cm"))
        
        # barplots of event counts
        counts <- as.numeric(metadata(x)$n_events)
        freqs <- round(counts/sum(counts)*100, 2)
        freq_labs <- paste0(counts, " (", freqs, "%)")
        freq_bars <- rowAnnotation(width=unit(2.5, "cm"), "n_events"=
                row_anno_barplot(x=counts, border=FALSE, axis=TRUE,
                    gp=gpar(fill="grey50", col="white"), bar_with=.75))
        freq_anno <- rowAnnotation(
            text=row_anno_text(freq_labs), 
            width=max_text_width(freq_labs))

        # heatmap of medians across antigens and samples
        heat_cols <- colorRampPalette(palette)(100)
        hm <- function(cell_fun) { Heatmap(
            med_exprs, heat_cols, "expression", 
            cell_fun=cell_fun, cluster_rows=row_clustering,
            heatmap_legend_param=list(color_bar="continuous"),
            column_names_gp=gpar(fontsize=8), row_names_gp=gpar(fontsize=8)) 
        }
        if (anno) {
            hm <- hm(cell_fun=function(j,i,x,y,...){
                grid.text(gp=gpar(fontsize=8), 
                    sprintf("%.2f", med_exprs[i,j]),x,y)})
        } else {
            hm <- hm(NULL)
        }
        p <- row_anno + hm + freq_bars + freq_anno
        
        if (!is.null(out_path)) {
            n <- ncol(med_exprs)
            out_nm <- file.path(out_path, "clustering_heatmap.pdf")
            pdf(out_nm, width=n/2, height=6); draw(p); dev.off()
        } else {
            p
        }
    }
)
