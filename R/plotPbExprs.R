#' @rdname plotPbExprs
#' @title Pseudobulk-level boxplot
#' 
#' @description 
#' Boxplot of aggregated marker data by sample or cluster, optionally 
#' colored and faceted by non-numeric cell metadata variables of interest.
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
#' @param color_by,group_by,shape_by 
#'   character string specifying a non-numeric cell metadata variable 
#'   to color, group and shape by, respectively; valid values are 
#'   \code{names(colData(x))} and \code{names(cluster_codes(x))} 
#'   if \code{\link{cluster}} has been run.
#' @param size_by logical specifying whether to scale point sizes by
#'   the number of cells in a given sample or cluster-sample instance;
#'   ignored when \code{geom = "boxes"}.
#' @param geom character string specifying whether 
#'   to include only points, boxplots or both.
#' @param jitter logical specifying whether to use \code{position_jitterdodge}
#'   in \code{geom_point} when \code{geom != "boxes"}.
#' @param ncol integer scalar specifying number of facet columns.
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
#' sce <- cluster(sce, verbose = FALSE)
#' 
#' # plot median expressions by sample & condition
#' # ...split by marker
#' plotPbExprs(sce, 
#'   shape_by = "patient_id",
#'   features = sample(rownames(sce), 6))
#' 
#' # ...split by cluster
#' plotPbExprs(sce, facet_by = "cluster_id", k = "meta6")
#' 
#' # plot median type-marker expressions by sample & cluster
#' plotPbExprs(sce, feature = "type", k = "meta6", 
#'   facet_by = "antigen", group_by = "cluster_id", color_by = "sample_id",
#'   size_by = TRUE, geom = "points", jitter = FALSE, ncol = 5)
#'   
#' # plot median state-marker expressions 
#' # by sample & cluster, split by condition
#' plotPbExprs(sce, k = "meta6", facet_by = "antigen", 
#'   group_by = "cluster_id", color_by = "condition", ncol = 7)
#' 
#' @import ggplot2
#' @importFrom dplyr across all_of group_by left_join mutate row_number
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment assay
#' @export

plotPbExprs <- function(x, k = "meta20", features = "state", 
    assay = "exprs", fun = c("median", "mean", "sum"), 
    facet_by = c("antigen", "cluster_id"), color_by = "condition", 
    group_by = color_by, shape_by = NULL, size_by = FALSE, 
    geom = c("both", "points", "boxes"), jitter = TRUE, ncol = NULL) {
    
    # check validity of input arguments
    fun <- match.arg(fun)
    geom <- match.arg(geom)
    facet_by <- match.arg(facet_by)
    stopifnot(is.logical(jitter), length(jitter) == 1)
    if (!is.null(ncol)) 
        stopifnot(is.numeric(ncol), length(ncol) == 1, ncol %% 1 == 0)
    if (facet_by == "cluster_id") {
        .check_sce(x, TRUE)
        k <- .check_k(x, k)
    } else .check_sce(x)
    .check_assay(x, assay)
    .check_cd_factor(x, color_by)
    .check_cd_factor(x, group_by)
    .check_cd_factor(x, shape_by)
    shapes <- .get_shapes(x, shape_by)
    if (is.null(shapes)) shape_by <- NULL
    x <- x[.get_features(x, features), ]

    # aggregation
    if (any(c(facet_by, group_by) == "cluster_id")) {
        x$cluster_id <- cluster_ids(x, k)
        by <- c("cluster_id", "sample_id")
    } else by <- "sample_id"
    ms <- .agg(x, by, fun, assay)
    df <- melt(ms, varnames = c("antigen", by[length(by)]))
    if (length(by) == 2) names(df)[ncol(df)] <- "cluster_id"
    x_var <- ifelse(facet_by == "antigen", group_by, "antigen")
    if (!is.null(df$cluster_id))
        df$cluster_id <- factor(df$cluster_id, levels(x$cluster_id))
    
    # add metadata information
    i <- match(df$sample_id, x$sample_id)
    j <- setdiff(names(colData(x)), c(names(df), "cluster_id"))
    df <- cbind(df, colData(x)[i, j, drop = FALSE])
    
    # add cell counts per sample(-cluster)
    ncs <- table(as.list(colData(x)[by]))
    ncs <- rep(c(t(ncs)), each = nrow(x))
    if (size_by) {
        size_by <- "n_cells" 
        df$n_cells <- ncs
    } else size_by <- NULL
    df <- df[ncs > 0, , drop = FALSE]

    ggplot(df, aes_string(x_var, "value", col = color_by)) +
        facet_wrap(facet_by, ncol = ncol, scales = "free_y") +
        (if (geom != "boxes") geom_point(
            alpha = 0.8, position = (if (jitter) {
                position_jitterdodge(jitter.width = 0.2, jitter.height = 0)
            } else "identity"),
            aes_string(fill = color_by, size = size_by, shape = shape_by))) +
        (if (geom != "points")
            geom_boxplot(alpha = 0.4, width = 0.8, fill = NA,
                outlier.color = NA, show.legend = FALSE)) +
        scale_shape_manual(values = shapes) + 
        scale_size_continuous(range = c(0.5, 3)) +
        guides(fill = "none", size = guide_legend(order = 3),
            shape = guide_legend(order = 2, override.aes = list(size = 3)),
            col = guide_legend(order = 1, 
                override.aes = list(alpha = 1, size = 3))) +
        ylab(paste(fun, ifelse(assay == "exprs", "expression", assay))) + 
        theme_bw() + theme(
            legend.key.height  =  unit(0.8, "lines"),
            axis.text = element_text(color = "black"),
            strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey", linewidth = 0.2)) +
        if (length(unique(c(x_var, color_by, group_by))) == 1) {
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        } else {
            theme(axis.text.x = element_text(
                angle = 45, hjust = 1, vjust = 1))
        }
}

#' @export
#' @rdname plotPbExprs
plotMedExprs <- function(x, 
    k = "meta20", features = "state",
    facet_by = c("antigen", "cluster_id"), 
    group_by = "condition", shape_by = NULL) {
    
    .Deprecated(
        old = "plotMedExprs",
        new = "plotPbExprs")
    
    plotPbExprs(x, k, features, 
        assay = "exprs", fun = "median", 
        facet_by, group_by, shape_by = shape_by)
}
