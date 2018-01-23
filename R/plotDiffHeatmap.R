# ==============================================================================
# Heatmap of median marker expressions across samples
# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
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
#' @return
#'
#' @references 
#' Nowicka M, Krieg C, Weber LM et al.
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return an object of class \code{HeatmapList{\link{ComplexHeatmap}}}.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' # specify contrasts
#' K <- matrix(c(0, 1), nrow=1, byrow=TRUE, dimnames=list("BCRXLvsRef"))
#' re <- diffAbundance(re, 12, K)
#' plotDiffHeatmap(re, 12, K)
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @import pheatmap SummarizedExperiment
#' @importFrom dplyr funs group_by summarize_all
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats hclust
#' @export
# ==============================================================================

setMethod(f="plotDiffHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, k, K, th=2.5, FDR_cutoff=.05, 
        type=c("abundance", "expr"),
        palette=brewer.pal(n=11, name="RdYlGn"), out_path=NULL) {
        
        cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
        counts <- as.data.frame.matrix(table(cluster_ids, sample_ids(x)))
        freqs <- t(t(counts) / colSums(counts))
        
        # arcsine-square-root transformation
        freqs_t <- as.data.frame.matrix(asin(sqrt(freqs)))
        # get significant clusters & sort
        nm <- rownames(K)
        type <- match.arg(type)
        adj_p <- setNames(
            metadata(x)[[paste0("diff_", type)]][[nm]][, paste0("adjp_", nm)],
            metadata(x)[[paste0("diff_", type)]][[nm]]$cluster_id)
        s <- names(which(sort(adj_p) < FDR_cutoff))
        adj_p <- adj_p[s]
        # Z-score normalization
        freqs_n <- Z_norm(freqs_t[s, ])
        
        # row & column labels
        row_labs <- paste0(s, " (", sprintf("%.02e", adj_p), ")")
        col_labs <- colnames(freqs_n)
        # row $ column annotations
        md <- metadata(x)[[1]]
        o <- order(md$condition)
        col_anno <- data.frame(condition=factor(md$condition)[o])
        rownames(col_anno) <- col_labs
        
        breaks <- seq(-th, th, length.out=101)
        lgd_breaks <- (-th):th
        pheatmap(freqs_n, border_color="white",
            gaps_col=as.numeric(table(col_anno$condition)),
            color=colorRampPalette(palette)(100), 
            breaks=breaks, legend_breaks=lgd_breaks,
            labels_row=row_labs, labels_col=col_labs,
            cluster_cols=FALSE, cluster_rows=FALSE,
            annotation_col=col_anno)
        
        # order samples by condition
        
    }
)
