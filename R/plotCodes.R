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
#' @param out_path character string. If specified, 
#'   output will be generated in this location.
#' @param verbose logical. Should information on progress be reported?
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
#' @import ggplot2
#' @importFrom cowplot get_legend plot_grid
#' @importFrom grDevices png
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom Rtsne Rtsne
#' @importFrom stats prcomp
#' @export

plotCodes <- function(x, k="meta20", out_path=NULL, verbose=TRUE) {
    
    # validity check
    .check_validity_of_k(x, k)
    stopifnot(c("cluster_codes", "SOM_codes") %in% names(metadata(x)))
    stopifnot(is.logical(verbose), length(verbose) == 1)
    if (!is.null(out_path))
        stopifnot(is.character(out_path), dir.exists(out_path))
    
    codes <- metadata(x)$SOM_codes
    if (verbose) message("o running tSNE...")
    tsne <- Rtsne(codes, pca=FALSE)
    if (verbose) message("o running PCA...")
    pca <- prcomp(codes, center=TRUE, scale.=FALSE)
    
    df <- data.frame(
        tSNE1=tsne$Y[, 1], tSNE2=tsne$Y[, 2],
        PCA1=pca$x[, 1], PCA2=pca$x[, 2])
    # get cluster IDs & sizes
    df$cluster_id <- cluster_codes(x)[, k]
    df$counts <- as.numeric(table(cluster_ids(x)))
    
    p <- ggplot(df, aes_string(color="cluster_id", size="counts")) +
        theme_classic() + theme(
            aspect.ratio=1, legend.position="top",
            panel.grid.minor=element_blank(),
            panel.grid.major=element_line(color='lightgrey', size=.2), 
            axis.title=element_text(face='bold'),
            axis.text=element_text(color="black"))
    
    # expand palette if more than 30 clusters
    nk <- nlevels(df$cluster_id)
    if (nk > 30) {
        cols <- colorRampPalette(.cluster_cols)(nk)
    } else {
        cols <- .cluster_cols[seq_len(nk)]
    }
    names(cols) <- levels(df$cluster_id)
    
    tsne_plot <- p + geom_point(aes_string(x="tSNE1", y="tSNE2")) + 
        scale_color_manual(values=cols, guide=FALSE) +
        scale_size(guide=FALSE)
    
    if (k > 10) n_row <- 2 else n_row <- 1
    pca_plot <- p + geom_point(aes_string(x="PCA1", y="PCA2")) +
        guides(color=guide_legend(override.aes=list(size=3), 
            order=1, nrow=n_row)) + scale_color_manual(values=cols) 
    
    p <- plot_grid(
        get_legend(pca_plot), 
        ncol = 1, rel_heights = c(1, 5),
        plot_grid(nrow = 1,
            tsne_plot + theme(legend.position = "none"), 
            pca_plot + theme(legend.position = "none")))
        
    if (!is.null(out_path)) {
        fn <- file.path(out_path, "codes_tsne+pca.png")
        ggsave(fn, p, width=9, height=5)
    } else {
        p
    }
}
