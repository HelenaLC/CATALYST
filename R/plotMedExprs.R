#' @rdname plotMedExprs
#' @title Plot median expressions
#' 
#' @description 
#' Plots median marker expressions across samples
#' computed on arcsinh-transformed intensities.
#'
#' @param x a \code{\link{SingleCellExperiment}{SingleCellExperiment}}.
#' @param k character string. Specifies the clustering to use.
#'   If \code{facet = "antigen"}, this argument will be ignored.
#' @param features specifies which features to include.
#'   Either a character string specifying a subset of features,
#'   or NULL for all features. When \code{rowData(x)$marker_class} 
#'   is specified, can be one of "type", "state", or "none".
#' @param facet \code{"antigen"} or \code{"cluster_id"}. 
#'   Note that the latter requires having run \code{\link{cluster}}.
#' @param group_by character string specifying 
#'   a \code{colData(x)} column to group samples by.
#' @param shape_by character string specifying 
#'   a \code{colData(x)} column to shape samples by.
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
#' plotMedExprs(sce, shape_by = "patient_id")
#' 
#' # run clustering
#' sce <- cluster(sce)
#' 
#' # plot median expressions across clusters
#' plotMedExprs(sce, facet = "cluster_id", k = "meta8")
#' 
#' @import ggplot2
#' @importFrom dplyr bind_rows
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment assay
#' @export

plotMedExprs <- function(x, 
    k="meta20", features="state",
    facet=c("antigen", "cluster_id"), 
    group_by="condition", shape_by=NULL) {
    
    .check_sce(x)
    .check_cd_factor(x, group_by)
    .check_cd_factor(x, shape_by)
    facet <- match.arg(facet)

    if (!is.null(shape_by)) {
        shapes <- c(16, 17, 15, 3, 7, 8) # default shapes
        if ((n <- nlevels(x[[shape_by]])) > 6) {
            if (n > 18) {
                message("At most 17 shapes are currently supported but ",
                    n, " are required. Setting 'shape_by' to NULL.")
                shape_by <- NULL
            } else {
                new <- setdiff(c(seq_len(16) - 1, 18), shapes)
                shapes <- c(shapes, new[seq_len(n - 6)])
            }
        }
    } else {
        shapes <- NULL
    }
    
    if (facet == "antigen") {
        ms <- .agg(x, by = "sample_id")
        ms <- melt(ms, varnames = c("antigen", "sample_id"))
    } else {
        k <- .check_validity_of_k(x, k)
        x$cluster_id <- cluster_ids(x, k)
        # subset features to use
        y <- x[.get_features(x, features), ]
        ms <- .agg(y, by = c("cluster_id", "sample_id"))
        ms <- lapply(ms, melt, varnames = c("antigen", "sample_id"))
        ms <- bind_rows(ms, .id = "cluster_id")
    }
    # add metadata information
    m <- match(ms$sample_id, ei(x)$sample_id)
    ms <- data.frame(ms, ei(x)[m, ])
    
    style <- list(ylab("median expression"),
        guides(color=guide_legend(override.aes=list(alpha=1))),
        geom_point(alpha=.75, aes_string(fill=group_by, shape=shape_by),
            position=position_jitterdodge(jitter.width=.25, jitter.height=0)),
        scale_shape_manual(values = shapes),
        scale_fill_manual(values = rep(NA, length(levels(ms[, group_by])))),
        geom_boxplot(alpha=.5, width=.75, fill=NA, outlier.color=NA),
        theme_bw(), theme(
            axis.text=element_text(color="black"),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_line(color="lightgrey", size=.25),
            strip.background=element_rect(fill="grey90", color=NA)))
    
    if (facet == "antigen") {
        ggplot(ms, aes_string(x=group_by, y="value", col=group_by)) +
            facet_wrap(facets="antigen", scales="free_y") + style + 
            theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    } else {
        ggplot(ms, aes_string(x="antigen", y="value", col=group_by)) + 
            facet_wrap(facets="cluster_id", scales="free_y", ncol=2) + style +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
    }
}
