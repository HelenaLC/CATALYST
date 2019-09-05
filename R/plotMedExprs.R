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
#' @param facet \code{"antigen"} or \code{"cluster_id"}. 
#'   Note that the latter requires having run \code{\link{cluster}}.
#' @param group_by character string specifying 
#'   a \code{colData(x)} column to group samples by.
#' @param shape_by character string specifying 
#'   a \code{colData(x)} column to shape samples by.
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
#' @importFrom dplyr group_by_ summarize_all
#' @importFrom matrixStats rowMedians
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment assay
#' @export

plotMedExprs <- function(x, k="meta20", 
    facet=c("antigen", "cluster_id"), 
    group_by="condition", shape_by=NULL) {
    
    .check_sce(x)
    facet <- match.arg(facet)
    stopifnot(is.character(group_by), 
        group_by %in% colnames(colData(x)),
        is.factor(x[[group_by]]))

    if (!is.null(shape_by)) {
        stopifnot(is.character(shape_by),
            shape_by %in% colnames(colData(x)), 
            is.factor(ei(x)[, shape_by]))
        shapes <- c(16, 17, 15, 3, 7, 8) # default shapes
        n <- nlevels(ei(x)[, shape_by])
        if (n > 6) {
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
        cs_by_s <- split(seq_len(ncol(x)), sample_ids(x))
        ms <- vapply(cs_by_s, function(cs) 
            rowMedians(assay(x, "exprs")[, cs, drop = FALSE]), 
            numeric(nrow(x)))
        ms <- data.frame(antigen = rownames(x), ms, check.names = FALSE)
        ms <- melt(ms, id.vars="antigen", variable.name="sample_id")
        
    } else {
        .check_validity_of_k(x, k)
        es <- assay(x, "exprs")[state_markers(x), ]
        ms <- data.frame(t(es), 
            sample_id=sample_ids(x), 
            cluster_id=cluster_ids(x, k)) %>% 
            group_by_(~sample_id, ~cluster_id) %>% 
            summarize_all(funs(median))
        ms <- melt(ms, variable.name="antigen",
            id.vars=c("sample_id", "cluster_id"))
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
