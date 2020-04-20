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
#'   valid values are \code{names(colData(x))} other than "sample_id"
#'   and "cluster_id"; NULL will use the first factor available.
#' @param shape_by character string specifying a non-numeric
#'   cell metadata columnd to shape by; valid values are 
#'   \code{names(colData(x))} other than "sample_id" and "cluster_id".
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
#' @export

plotAbundances <- function(x, k = "meta20", 
    by = c("sample_id", "cluster_id"), 
    group_by = NULL, shape_by = NULL,
    k_pal = CATALYST:::.cluster_cols) {
    
    # check validity of input arguments
    by <- match.arg(by)
    .check_sce(x, TRUE)
    k <- .check_k(x, k)
    .check_cd_factor(x, group_by)
    .check_cd_factor(x, shape_by)
    .check_pal(k_pal)
    
    shapes <- .get_shapes(x, shape_by)
    if (is.null(shapes)) shape_by <- NULL
    
    valid <- setdiff(names(colData(x)), c("sample_id", "cluster_id"))
    if (length(valid) == 0)
        stop("No factors to group by. Metadata should contain", 
            " at least one column other than 'file' and 'id'.")
    if (is.null(group_by)) group_by <- valid[1]
    
    # ramp cluster color palette
    if (by == "sample_id") {
        nk <- nlevels(cluster_ids(x, k))
        if (length(k_pal) < nk)
            k_pal <- colorRampPalette(k_pal)(nk)
    }
    
    # get frequencies by cluster & sample
    ns <- table(cluster_ids(x, k), sample_ids(x))
    fq <- prop.table(ns, 2) * 100
    df <- melt(fq, varnames = c("cluster_id", "sample_id"))
    
    # add relevant cell metadata
    m <- match(df$sample_id, x$sample_id)
    for (i in c(shape_by, group_by))
        df[[i]] <- x[[i]][m]
    
    # specify shared aesthetics
    p <- ggplot(df, aes_string(y = "value")) +
        labs(x = NULL, y = "Proportion [%]") + 
        theme_bw() + theme(
            panel.grid = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.key.height  =  unit(0.8, "lines"))
    
    switch(by,
        sample_id = p + 
            facet_wrap(group_by, scales = "free_x") +
            geom_bar(
                aes_string(x = "sample_id", fill = "factor(cluster_id)"), 
                position = "fill", stat = "identity") +
            scale_fill_manual("cluster_id", values = k_pal) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
            theme(
                panel.border = element_blank(),
                panel.spacing.x = unit(1, "lines"))
        ,
        cluster_id = p + 
            facet_wrap("cluster_id", scales = "free_y", ncol = 4) +
            geom_boxplot(
                aes_string(x = group_by, color = group_by, fill = group_by),
                position = position_dodge(), alpha = 0.2, outlier.color = NA) + 
            geom_point(
                aes_string(x = group_by, col = group_by, shape = shape_by),
                position = position_jitter(width = 0.2)) +
            scale_shape_manual(values = shapes) + guides(fill = FALSE, 
                shape = guide_legend(override.aes = list(size = 3))) 
    )
}
