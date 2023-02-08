#' @rdname plotClusterExprs
#' @title Plot expression distributions by cluster
#' 
#' @description Plots smoothed densities of marker intensities by cluster.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k character string specifying which clustering to use;
#'   valid values are \code{names(cluster_codes(x))}.
#' @param features a character vector specifying 
#'   which antigens to include; valid values are
#'   \code{"type"/"state"} for \code{type/state_markers(x)} 
#'   if \code{rowData(x)$marker_class} have been specified; 
#'   a subset of \code{rownames(x)}; NULL to use all features.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Crowell HL, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return a \code{\link{ggplot}} object.
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' plotClusterExprs(sce, k = "meta8")
#' 
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay colData
#' @export

plotClusterExprs <- function(x, 
    k = "meta20", features = "type") {
    
    # check validity of input arguments
    .check_sce(x, TRUE)
    k <- .check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
    features <- .get_features(x, features)

    # order clusters according to hierarchical 
    # clustering on median feature expressions 
    ms <- t(.agg(x[features, ], "cluster_id", "median"))
    d <- dist(ms, method="euclidean")
    o <- hclust(d, method="average")$order
    
    # construct data.frame of expression matrix include cell metadata
    cd <- colData(x)
    es <- assay(x[features, ], "exprs")
    df <- data.frame(t(es), cd, check.names = FALSE)
    df <- melt(df, 
        id.vars = names(cd),
        variable.name = "antigen", 
        value.name = "expression")
    # add average across all clusters as reference
    df$avg <- "no"
    avg <- df
    avg$cluster_id <- "avg"
    avg$avg <- "yes"
    df <- rbind(df, avg)
    
    # compute cluster frequencies
    fq <- tabulate(x$cluster_id) / ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    
    # reorder clusters
    df$cluster_id <- factor(df$cluster_id, 
        levels = rev(c("avg", levels(x$cluster_id)[o])),
        labels = rev(c("average", paste0(names(fq), " (", fq, "%)")[o])))
    
    ggplot(df, aes_string(
        x = "expression", y = "cluster_id", 
        col = "avg", fill = "avg")) + 
        facet_wrap(~antigen, scales = "free_x", nrow = 2) + 
        geom_density_ridges(alpha = 0.2) + 
        theme_ridges() + theme(
            legend.position = "none",
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"))
}