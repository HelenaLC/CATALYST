# ==============================================================================
# Barplot of relative population abundances across samples & clusters
# ------------------------------------------------------------------------------
#' @rdname plotAbundances
#' @title Population frequencies across samples & clusters
#' 
#' @description 
#' Plots the relative population abundances of the specified clustering.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k character string specifying which clustering to use;
#'   valid values are \code{names(cluster_codes(x))}.
#' @param by a character string specifying whether 
#'   to plot frequencies by samples or clusters.
#' @param group_by character string specifying a non-numeric
#'   cell metadata columnd to group by (determines the color coding);
#'   valid values are \code{names(colData(x))} 
#'   other than "sample_id" and "cluster_id".
#' @param shape_by character string specifying a non-numeric
#'   cell metadata columnd to shape by; valid values are 
#'   \code{names(colData(x))} other than "sample_id" and "cluster_id".
#' @param col_clust for \code{by = "sample_id"}, 
#'   specifies whether to hierarchically cluster samples 
#'   and reorder them accordingly. When \code{col_clust = FALSE}, 
#'   samples are ordered according to \code{levels(x$sample_id)} 
#'   (or alphabetically, when \code{x$sample_id} is not a factor).
#' @param distance character string specifying the distance metric 
#'  to use for sample clustering; passed to \code{\link[stats]{dist}} 
#' @param linkage character string specifying the agglomeration method 
#'  to use for sample clustering; passed to \code{\link[stats]{hclust}}.
#' @param k_pal character string specifying the cluster 
#'   color palette; ignored when \code{by = "cluster_id"}. 
#'   If less than \code{nlevels(cluster_ids(x, k))} 
#'   are supplied, colors will be interpolated via 
#'   \code{\link[grDevices:colorRamp]{colorRampPalette}}.
#'   
#' @return a \code{ggplot} object.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Crowell HL, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' # plot relative population abundances 
#' # by sample & cluster, respectively
#' plotAbundances(sce, k = "meta12")                  
#' plotAbundances(sce, k = "meta8", by = "cluster_id") 
#' 
#' # use custom cluster color palette
#' plotAbundances(sce, k = "meta10", 
#'   k_pal = c("lightgrey", "cornflowerblue", "navy"))
#' 
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @export

plotAbundances <- function(x, k = "meta20", 
    by = c("sample_id", "cluster_id"), 
    group_by = "condition", shape_by = NULL,
    col_clust = TRUE, 
    distance = c(
        "euclidean", "maximum", "manhattan", 
        "canberra", "binary", "minkowski"), 
    linkage = c(
        "average", "ward.D", "single", "complete", 
        "mcquitty", "median", "centroid", "ward.D2"),
    k_pal = .cluster_cols) {
    
    # check validity of input arguments
    by <- match.arg(by)
    .check_sce(x, TRUE)
    k <- .check_k(x, k)
    .check_cd_factor(x, group_by)
    .check_cd_factor(x, shape_by)
    .check_pal(k_pal)
    linkage <- match.arg(linkage)
    distance <- match.arg(distance)
    stopifnot(is.logical(col_clust), length(col_clust) == 1)
    
    shapes <- .get_shapes(x, shape_by)
    if (is.null(shapes)) shape_by <- NULL
    
    # ramp cluster color palette
    if (by == "sample_id") {
        nk <- nlevels(cluster_ids(x, k))
        if (length(k_pal) < nk)
            k_pal <- colorRampPalette(k_pal)(nk)
    }
    
    # get frequencies by cluster & sample
    ns <- table(
        cluster_id = cluster_ids(x, k), 
        sample_id = sample_ids(x))
    fq <- prop.table(ns, 2) * 100
    df <- as.data.frame(fq)
    
    # add relevant cell metadata
    m <- match(df$sample_id, x$sample_id)
    for (i in c(shape_by, group_by))
        df[[i]] <- x[[i]][m]
    
    if (by == "sample_id" && col_clust 
        && length(unique(df$sample_id)) > 1) {
        d <- dist(t(fq), distance)
        h <- hclust(d, linkage)
        o <- colnames(fq)[h$order]
        df$sample_id <- factor(df$sample_id, o)
    }

    # specify shared aesthetics
    p <- ggplot(df, aes_string(y = "Freq")) +
        labs(x = NULL, y = "Proportion [%]") + 
        theme_bw() + theme(
            panel.grid = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.key.height  =  unit(0.8, "lines"))
    
    switch(by,
        sample_id = p + (if (!is.null(group_by)) 
            facet_wrap(group_by, scales = "free_x")) +
            geom_bar(
                aes_string(x = "sample_id", fill = "cluster_id"), 
                position = "fill", stat = "identity") +
            scale_fill_manual("cluster_id", values = k_pal) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
            theme(
                panel.border = element_blank(),
                panel.spacing.x = unit(1, "lines"))
        ,
        cluster_id = {
            p <- p + scale_shape_manual(values = shapes) + guides(
                col = guide_legend(order = 1, override.aes = list(size = 3)),
                shape = guide_legend(override.aes = list(size = 3)))
            if (is.null(group_by)) {
                p + geom_boxplot(aes_string(x = "cluster_id"), alpha = 0.2,
                    position = position_dodge(), outlier.color = NA) + 
                    geom_point(aes_string("cluster_id", shape = shape_by),
                        position = position_jitter(width = 0.2))
            } else {
                p + facet_wrap("cluster_id", scales = "free_y", ncol = 4) +
                    geom_boxplot(aes_string(x = group_by, 
                        color = group_by, fill = group_by),
                        position = position_dodge(), alpha = 0.2, 
                        outlier.color = NA, show.legend = FALSE) + 
                    geom_point(aes_string(x = group_by, 
                        col = group_by, shape = shape_by),
                        position = position_jitter(width = 0.2))
            }
        }
    )
}
