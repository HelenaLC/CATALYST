#' @rdname plotCodes
#' @title tSNE and PCA on SOM codes
#' 
#' @description Plots the tSNE and PCA representing the SOM codes as inferred
#' by \pkg{FlowSOM}. Sizes are scaled to the total number of events assigned 
#' to each cluster, and points are color according to their cluster ID upon 
#' \pkg{ConsensusClusterPlus} metaclustering into \code{k} clusters.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param k 
#'   character string. 
#'   Specifies the clustering to use for color coding.
#' @param out_path 
#'   character string. If specified, output will be generated in this location.
#' @param verbose 
#'   logical. Specifies whether information on progress should be reported.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
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
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run clustering
#' re <- cluster(re)
#' plotCodes(re)
#' 
#' @import ggplot2 Rtsne
#' @import Rtsne  
#' @importFrom dplyr funs group_by summarize_all
#' @importFrom grDevices png
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom stats prcomp
# ------------------------------------------------------------------------------

setMethod(f="plotCodes", 
    signature=signature(x="daFrame"), 
    definition=function(x, k="meta20", out_path=NULL, verbose=TRUE) {
        
        # validity check
        check_validity_of_k(x, k)
        
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
        n_clusters <- nlevels(df$cluster_id)
        if (n_clusters > 30) {
            cols <- colorRampPalette(cluster_cols)(n_clusters)
        } else {
            cols <- cluster_cols[seq_len(n_clusters)]
        }
        names(cols) <- levels(df$cluster_id)
        
        tsne_plot <- p + geom_point(aes_string(x="tSNE1", y="tSNE2")) + 
            scale_color_manual(values=cols, guide=FALSE) +
            scale_size(guide=FALSE)
        
        if (k > 10) n_row <- 2 else n_row <- 1
        pca_plot <- p + geom_point(aes_string(x="PCA1", y="PCA2")) +
            guides(color=guide_legend(override.aes=list(size=3), 
                nrow=n_row)) + scale_color_manual(values=cols) 
        
        # store legend
        legend <- get_legend(pca_plot)
        legend_height <- sum(legend$heights)
        pca_plot <- pca_plot + theme(legend.position="none")
        
        # align heights
        gA <- ggplotGrob(tsne_plot)
        gB <- ggplotGrob(pca_plot)
        gA$widths <- gB$widths
        
        g <- arrangeGrob(legend, tsne_plot, pca_plot,
            layout_matrix=matrix(c(1,2,1,3), nrow=2),
            heights=unit.c(legend_height, unit(1, "npc")-legend_height))
        if (!is.null(out_path)) {
            ggsave(file.path(out_path, "codes_tsne+pca.png"), 
                g, width=9, height=5)
        } else {
            grid.arrange(g)
        }
    }
)
