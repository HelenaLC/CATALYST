#' @rdname plotAggExprs
#' @title Boxplot of aggregated marker data
#' 
#' @description 
#' Boxplots of aggregated marker data by sample or cluster,
#' colored by non-numeric cell metadata variable of interest.
#'
#' @param x a \code{\link{SingleCellExperiment}{SingleCellExperiment}}.
#' @param k character string specifying which clustering to use;
#'   values values are \code{names(cluster_codes(x))}.
#'   Ignored if \code{facet_by = "antigen"}.
#' @param features character vector specifying 
#'   which features to include; valid values are 
#'   \code{"type"/"state"} for \code{type/state_markers(x)} 
#'   if \code{rowData(x)$marker_class} have been specified; 
#'   a subset of \code{rownames(x)}; NULL to use all features.
#' @param assay character string specifying which assay data 
#'   to use; valid values are \code{assayNames(x)}.
#' @param fun character string specifying the summary statistic to use.
#' @param facet_by \code{"antigen"} or \code{"cluster_id"}; 
#'   the latter requires having run \code{\link{cluster}}.
#' @param group_by character string specifying a non-numeric cell metadata 
#'   variable to group samples by; valid values are \code{names(colData(x))}.
#' @param shape_by character string specifying a non-numeric cell metadata 
#'   variable to shape samples by; valid values are \code{names(colData(x))}.
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
#' # construct SCE
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # plot median expressions
#' plotAggExprs(sce, 
#'   shape_by = "patient_id",
#'   features = sample(rownames(sce), 6))
#' 
#' # run clustering
#' sce <- cluster(sce, verbose = FALSE)
#' 
#' # plot median expressions across clusters
#' p <- plotAggExprs(sce, facet_by = "cluster_id", k = "meta8")
#' p$facet$params$ncol <- 4; p
#' 
#' # change colors & facetting layout
#' library(ggplot2)
#' p$facet$params$ncol <- 4
#' p + scale_color_manual(values = c(
#'   Ref = "royalblue", BCRXL = "orange"))
#' 
#' @import ggplot2
#' @importFrom dplyr bind_rows
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment assay
#' @export

plotAggExprs <- function(x, 
    k = "meta20", features = "state",
    assay = "exprs", fun = c("median", "mean", "sum"),
    facet_by = c("antigen", "cluster_id"), 
    group_by = "condition", shape_by = NULL) {
    
    # check validity of input arguments
    fun <- match.arg(fun)
    facet_by <- match.arg(facet_by)
    if (facet_by == "antigen") {
        .check_sce(x)
    } else {
        .check_sce(x, TRUE)
        k <- .check_k(x, k)
    }
    .check_assay(x, assay)
    .check_cd_factor(x, group_by)
    .check_cd_factor(x, shape_by)
    features <- .get_features(x, features)
    shapes <- .get_shapes(x, shape_by)
    if (is.null(shapes)) shape_by <- NULL
    
    x <- x[features, ]
    if (facet_by == "antigen") {
        ms <- .agg(x, "sample_id", fun, assay)
        df <- melt(ms, varnames = c("antigen", "sample_id"))
        # aesthetics
        x_var <- group_by
        thm <-  theme(
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    } else {
        x$cluster_id <- cluster_ids(x, k)
        ms <- .agg(x, c("cluster_id", "sample_id"), fun, assay)
        dfs <- lapply(ms, melt, varnames = c("antigen", "sample_id"))
        df <- bind_rows(dfs, .id = "cluster_id")
        df$cluster_id <- factor(df$cluster_id, 
            levels = levels(x$cluster_id))
        # aesthetics
        x_var <- "antigen"
        thm <-  theme(axis.text.x = element_text(
            angle = 45, hjust = 1, vjust = 1))
    }
    # add metadata information
    m <- match(df$sample_id, x$sample_id)
    for (i in c(group_by, shape_by))
        df[[i]] <- x[[i]][m]
    
    ylab <- paste(fun, ifelse(assay == "exprs", "expression", assay))
    ggplot(df, aes_string(x_var, "value", col = group_by)) +
        facet_wrap(facet_by, scales = "free_y") +
        geom_point(alpha = 0.8, 
            aes_string(fill = group_by, shape = shape_by),
            position = position_jitterdodge(
                jitter.width = 0.2, jitter.height = 0)) +
        geom_boxplot(alpha = 0.4, width = 0.8, 
            fill = NA, outlier.color = NA, show.legend = FALSE) +
        scale_shape_manual(values = shapes) + guides(
            shape = guide_legend(override.aes = list(size = 3)),
            col = guide_legend(override.aes = list(alpha = 1, size = 3))) +
        ylab(ylab) + theme_bw() + thm + theme(
            legend.key.height  =  unit(0.8, "lines"),
            axis.text = element_text(color = "black"),
            strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey", size = 0.2))
}

#' @export
#' @rdname plotAggExprs
plotMedExprs <- function(x, 
    k = "meta20", features = "state",
    facet_by = c("antigen", "cluster_id"), 
    group_by = "condition", shape_by = NULL) {
    
    .Deprecated(
        old = "plotMedExprs",
        new = "plotAggExprs")
    
    plotAggExprs(x, k, features, 
        assay = "exprs", fun = "median", 
        facet_by, group_by, shape_by)
}
