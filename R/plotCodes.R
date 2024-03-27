#' @rdname plotCodes
#' @title tSNE and PCA on SOM codes
#' 
#' @description Plots the tSNE and PCA representing the SOM codes as inferred
#' by \pkg{FlowSOM}. Sizes are scaled to the total number of events assigned 
#' to each cluster, and points are color according to their cluster ID upon 
#' \pkg{ConsensusClusterPlus} metaclustering into \code{k} clusters.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k character string. Specifies the clustering to use for color coding.
#' @param k_pal character string specifying the cluster color palette;
#'   If less than \code{nlevels(cluster_ids(x, k))} are supplied, colors will 
#'   be interpolated via \code{\link[grDevices:colorRamp]{colorRampPalette}}.
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
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' plotCodes(sce, k = "meta14")
#' 
#' # use custom cluster color palette
#' plotCodes(sce, k = "meta12",
#'   k_pal = c("lightgrey", "cornflowerblue", "navy"))
#' 
#' @import ggplot2
#' @importFrom cowplot get_legend plot_grid
#' @importFrom grDevices colorRampPalette
#' @importFrom Rtsne Rtsne
#' @importFrom stats prcomp
#' @importFrom S4Vectors metadata
#' @export

plotCodes <- function(x, k = "meta20", k_pal = .cluster_cols) {
    
    # check validity of input arguments
    .check_sce(x, TRUE)
    .check_pal(k_pal)
    k <- .check_k(x, k)
    
    # ramp cluster color palette
    nk <- nlevels(cluster_ids(x, k))
    if (length(k_pal) < nk)
        k_pal <- colorRampPalette(k_pal)(nk)
    
    # run t-SNE & PCA on SOM codes
    codes <- metadata(x)$SOM_codes
    tsne <- Rtsne(codes, pca = FALSE)
    pca <- prcomp(codes, center = TRUE, scale. = FALSE)

    # construct data.frame of t-SNE & PCA coordinates
    df <- data.frame(tsne$Y, pca$x[, c(1, 2)]) 
    colnames(df) <- c("tsne1", "tsne2", "pc1", "pc2")
    
    # add cluster IDs & sizes
    df$cluster_id <- cluster_codes(x)[, k]
    df$counts <- as.numeric(table(cluster_ids(x)))
    
    # specify shared aesthetics
    p <- ggplot(df, aes_string(col = "cluster_id", size = "counts")) +
        theme_classic() + theme(
            aspect.ratio = 1, 
            legend.position = "top",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey", linewidth = 0.2), 
            axis.text=element_text(color = "black"))

    ps <- list(
        p + geom_point(aes_string("tsne1", "tsne2")) + 
            labs(x = "t-SNE dim. 1", y = "t-SNE dim. 2") +
            scale_color_manual(values = k_pal, guide = "none") +
            scale_size(guide = "none")
        ,
        p + geom_point(aes_string("pc1", "pc2")) +
            labs(x = "1st PC", y = "2nd PC") +
            scale_color_manual(values = k_pal) +
            guides(col = guide_legend(
                override.aes = list(size = 3),
                order = 1, nrow = ifelse(nk > 10, 2, 1)))
    )
    
    # arrange plots side-by-side
    suppressWarnings(lgd <- get_legend(ps[[2]]))
    ps <- lapply(ps, "+", theme(legend.position = "none"))
    plot_grid(
        lgd, ncol = 1, rel_heights = c(1, 5),
        plot_grid(plotlist = ps, nrow = 1))
}

