#' @rdname pbMDS
#' @title Pseudobulk-level MDS plot
#' 
#' @description Pseudobulk-level Multi-Dimensional Scaling (MDS) 
#' plot computed on median marker expressions in each sample.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param by character string specifying whether to aggregate 
#'   by \code{sample_id}, \code{cluster_id} or \code{both}.
#' @param k character string specifying which clustering to use when 
#'   \code{by != "sample_id"}; valid values are \code{names(cluster_codes(x))}.
#' @param dims two numeric scalars indicating which dimensions to plot.
#' @param features character string specifying which features to include
#'   for computation of reduced dimensions; valid values are 
#'   \code{"type"/"state"} for \code{type/state_markers(x)} 
#'   if \code{rowData(x)$marker_class} have been specified; 
#'   a subset of \code{rownames(x)}; NULL to use all features.
#' @param assay character string specifying which assay data to use;
#'   valid values are \code{assayNames(x)}.
#' @param fun character string specifying which summary statistic to use.
#' @param color_by character string specifying a 
#'   non-numeric cell metadata column to color by; 
#'   valid values are \code{names(colData(x))}.
#' @param label_by character string specifying a 
#'   non-numeric cell metadata column to label by; 
#'   valid values are \code{names(colData(x))}.
#' @param shape_by character string specifying a 
#'   non-numeric cell metadata column to shape by; 
#'   valid values are \code{names(colData(x))}.
#' @param size_by logical specifying whether points should be 
#'   sized by the number of cells that went into aggregation; i.e., 
#'   the size of a give sample, cluster or cluster-sample instance.
#' @param pal character vector of colors to use; 
#'   NULL for default \code{ggplot2} colors.
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
#' sce <- cluster(sce)
#' 
#' # sample-level pseudobulks
#' # including state-markers only
#' pbMDS(sce, by = "sample_id", features = "state")
#' 
#' # cluster-level pseudobulks
#' # including type-features only
#' pbMDS(sce, by = "cluster_id", features = "type")
#' 
#' # pseudobulks by cluster-sample 
#' # including all features
#' pbMDS(sce, by = "both", k = "meta12", 
#'   shape_by = "condition", size_by = TRUE)
#' 
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom scater calculateMDS
#' @importFrom SummarizedExperiment colData
#' @export

pbMDS <- function(x,
    by = c("sample_id", "cluster_id", "both"), k = "meta20", dims = c(1, 2),
    features = NULL, assay = "exprs", fun = c("median", "mean", "sum"), 
    color_by = switch(by, sample_id = "condition", "cluster_id"),
    label_by = if (by == "sample_id") "sample_id" else NULL, 
    shape_by = NULL, size_by = is.null(shape_by),
    pal = if (color_by == "cluster_id") .cluster_cols else NULL) {
    
    # check validity of input arguments
    by <- match.arg(by)
    fun <- match.arg(fun)
    args <- as.list(environment())
    .check_args_pbMDS(args)
    
    if (by != "sample_id") 
        x$cluster_id <- cluster_ids(x, k)
    by <- switch(by, both = c("cluster_id", "sample_id"), by)
    
    # aggregate & run MDS
    x <- x[.get_features(x, features), ]
    pbs <- .agg(x, by, fun, assay)
    if (is.list(pbs))
        pbs <- do.call("cbind", pbs)
    mds <- calculateMDS(pbs, ncomponents = max(dims))
    
    # construct data.frame for plotting
    df <- data.frame(mds[, dims])
    colnames(df) <- c("x", "y")
    if (length(by) == 1) {
        df[[by]] <- factor(colnames(pbs), levels(x[[by]]))
    } else {
        ns <- length(sids <- levels(x$sample_id))
        nk <- length(kids <- levels(x$cluster_id))
        df$sample_id <- factor(rep(sids, nk), sids)
        df$cluster_id <-  factor(rep(kids, each = ns), kids)
    }
    
    # add sample metadata
    if (!isTRUE(by == "cluster_id")) {
        m <- match(df$sample_id, x$sample_id)
        i <- setdiff(names(colData(x)), names(df))
        df <- cbind(df, colData(x)[m, i, drop = FALSE])
    }
    
    # add instance cell counts
    if (size_by) {
        size_by <- "n_cells" 
        df$n_cells <- c(t(table(as.list(colData(x)[by]))))
    } else size_by <- NULL
    
    # get number of legend columns to use
    ncol <- ifelse(!is.null(color_by) && nlevels(df[[color_by]]) > 10, 2, 1)
    
    ggplot(df, aes(.data$x, .data$y, 
        col = if (!is.null(color_by)) .data[[color_by]], 
        shape = if (!is.null(shape_by)) .data[[shape_by]])) + 
        geom_point(alpha = 0.8, 
            aes(size = if (!is.null(size_by)) .data[[size_by]])) + 
        (if (!is.null(label_by)) geom_label_repel(
            aes(label = .data[[label_by]]), show.legend = FALSE)) + 
        (if (!is.null(pal)) scale_color_manual(values = pal)) +
        scale_shape_manual(values = .get_shapes(x, shape_by)) +
        guides(
            col = guide_legend(order = 1, ncol = ncol, 
                override.aes = list(alpha = 1, size = 3)),
            shape = guide_legend(order = 2, override.aes = list(size = 3)),
            size = guide_legend(order = 3)) + 
        labs(
            x = paste("MDS dim.", dims[1]), 
            y = paste("MDS dim.", dims[2]),
            col=color_by, shape=shape_by, size=size_by) +
        coord_equal() + theme_linedraw() + theme(
            panel.grid.minor = element_blank(),
            legend.key.height  =  unit(0.8, "lines"))
}