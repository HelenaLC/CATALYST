#' @rdname plotDR
#' @title Plot reduced dimensions
#' 
#' @description Dimension reduction plot colored 
#' by expression, cluster, sample or group ID.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param dr character string specifying which dimension reduction to use. 
#'   Should be one of \code{reducedDimNames(x)}; default to the 1st available.
#' @param color_by character string corresponding to a
#'   \code{colData(x)} column. Specifies the color coding.
#' 
#' @author Helena Lucia Crowell
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' library(scater)
#' sce <- runUMAP(sce)
#' 
#' # color by pS6 expression, split by group
#' plotDR(sce, color_by = "pS6") +
#'   facet_wrap("condition")
#' 
#' # color by 8 metaclusters, split by sample
#' plotDR(sce, color_by = "meta8") + 
#'   facet_wrap("sample_id", ncol = 4)
#' 
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom scales hue_pal
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom SummarizedExperiment assay colData
#' @export

plotDR <- function(x, dr = NULL, color_by = "condition") {
    
    # check validity of input arguments
    kids <- names(cluster_codes(x))
    u <- c(rownames(x), colnames(colData(x)), kids)
    stopifnot(is.character(color_by), length(color_by) == 1, color_by %in% u)
    .check_sce(x, color_by %in% kids)
    stopifnot(length(reducedDimNames(x)) != 0)
    if (is.null(dr)) dr <- reducedDimNames(x)[1]
    stopifnot(is.character(dr), length(dr) == 1, dr %in% reducedDimNames(x))
    
    # construct data.frame
    xy <- reducedDim(x, dr)
    df <- data.frame(colData(x), xy)
    if (color_by %in% rownames(x))
        df[[color_by]] <- assay(x, "exprs")[color_by, ]
    if (color_by %in% kids)
        df[[color_by]] <- cluster_ids(x, color_by)
    
    if (is(df[[color_by]], "numeric")) {
        col <- scale_color_viridis_c()
    } else {
        if (color_by %in% kids) {
            col <- .cluster_cols
            if (n <- nlevels(df[[color_by]]) > length(col))
                col <- colorRampPalette(col)(n)
        } else {
            col <- hue_pal()(nlevels(df[[color_by]]))
        }
        col <- list(
            scale_color_manual(values = col),
            theme(legend.key.size = unit(2, "mm")),
            guides(color = guide_legend(
                override.aes = list(alpha = 1, size = 2))))
    }
    
    ggplot(df, aes_string(x = "X1", y = "X2", col = color_by)) +
        geom_point(size = 0.4, alpha = 0.8) + col + 
        labs(x = paste(dr, "dim. 1"), y = paste(dr, "dim. 2")) +
        theme_minimal() + theme(aspect.ratio = 1,
            panel.grid.minor = element_blank(),
            axis.text = element_text(color = "black"))
}