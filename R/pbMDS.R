#' @rdname pbMDS
#' @title Pseudobulk-level MDS plot
#' 
#' @description Pseudobulk-level Multi-Dimensional Scaling (MDS) 
#' plot computed on median marker expressions in each sample.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param color_by character string specifying a 
#'   non-numeric cell metadata column to color by; 
#'   valid values are \code{names(colData(x))}.
#' @param label_by character string specifying a 
#'   non-numeric cell metadata column to label by; 
#'   valid values are \code{names(colData(x))}.
#' @param shape_by character string specifying a 
#'   non-numeric cell metadata column to shape by; 
#'   valid values are \code{names(colData(x))}.
#' @param assay character string specifying which assay data to use;
#'   valid values are \code{assayNames(x)}.
#' @param fun character string specifying 
#'   which function to use as summary statistic.
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
#' pbMDS(sce)
#' 
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom limma plotMDS
#' @importFrom matrixStats rowMedians
#' @importFrom methods is
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay assayNames colData
#' @export

pbMDS <- function(x, 
    color_by = "condition", label_by = "sample_id", shape_by = NULL,
    assay = "exprs", fun = c("median", "mean", "sum")) {
    # check validity of input arguments
    .check_sce(x)
    .check_assay(x, assay)
    .check_cd_factor(x, color_by)
    .check_cd_factor(x, label_by)
    .check_cd_factor(x, shape_by)
    fun <- match.arg(fun)
    
    # compute medians across samples
    y <- .agg(x, "sample_id", fun, assay)
    y <- y[, colSums(y) != 0]
    
    # get MDS coordinates
    mds <- plotMDS(y, plot = FALSE)
    df <- data.frame(mds1 = mds$x, mds2 = mds$y)
    
    # add relevant cell metadata
    m <- match(rownames(df), x$sample_id)
    for (i in c(color_by, label_by, shape_by))
        df[[i]] <- x[[i]][m]
    
    ggplot(df, aes_string("mds1", "mds2", col = color_by, shape = shape_by)) + 
        geom_point(alpha = 0.8, size = 2) + 
        geom_label_repel(aes_string(label = label_by), show.legend = FALSE) + 
        guides(shape = guide_legend(override.aes = list(size = 3)),
            col = guide_legend(override.aes = list(alpha = 1, size = 3))) +
        labs(x = "MDS dim. 1", y = "MDS dim. 2") +
        coord_equal() + theme_minimal() + theme(
            legend.key.height  =  unit(0.8, "lines"),
            axis.text = element_text(color = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey", size = 0.2))
}
