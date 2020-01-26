#' @rdname plotMDS
#' @title MDS plot
#' 
#' @description 
#' Multi-dimensional scaling (MDS) plot on median marker expressions.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param color_by character string corresponding to a
#'   \code{colData(x)} column. Specifies the color coding.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Crowell HL, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' plotMDS(sce)
#' 
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom limma plotMDS
#' @importFrom matrixStats rowMedians
#' @importFrom methods is
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay colData
#' @export

plotMDS <- function(x, color_by="condition") {
    .check_sce(x)
    stopifnot(is.character(color_by), length(color_by) == 1, 
        color_by %in% names(colData(x)))
    
    # compute medians across samples
    cs_by_s <- split(seq_len(ncol(x)), droplevels(x$sample_id))
    es <- as.matrix(assay(x, "exprs"))
    ms <- vapply(cs_by_s, function(cs)
        rowMedians(es[, cs, drop = FALSE]),
        numeric(nrow(x)))
    rownames(ms) <- rownames(x)
    
    # get MDS coordinates
    mds <- limma::plotMDS(ms, plot=FALSE)
    df <- data.frame(MDS1=mds$x, MDS2=mds$y)
    
    # add metadata information
    md <- metadata(x)$experiment_info
    m <- match(rownames(df), md$sample_id)
    df <- data.frame(df, md[m, ])
    
    ggplot(df, aes_string(x="MDS1", y="MDS2", col=color_by)) + 
        geom_label_repel(aes_string(label="sample_id"), 
            show.legend=FALSE) + geom_point(alpha=.8, size=1.2) + 
        guides(col=guide_legend(overide.aes=list(alpha=1, size=3))) +
        theme_void() + theme(aspect.ratio=1,
            panel.grid.minor=element_blank(),
            panel.grid.major=element_line(color='lightgrey', size=.25), 
            axis.title=element_text(face='bold'),
            axis.title.y=element_text(angle=90),
            axis.text=element_text())
}