#' @rdname plotClusterExprs
#' @title Plot expression distributions by cluster
#' 
#' @description 
#' Plots smoothed densities of arcsinh-transformed 
#' marker intensities by cluster.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k 
#'   character string. Specifies the clustering to use.
#' @param features
#'   character string specifying which features to include. Defaults to NULL 
#'   (= all features). Alternatively, if the \code{colData(x)$marker_class} 
#'   column is specified, can be one of "type", "state", or "none".
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
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
#' plotClusterExprs(sce)
#' 
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay colData
#' @export

plotClusterExprs <- function(x, k="meta20", features=NULL) {
    
    # check validity of input arguments
    .check_sce(x, TRUE)
    k <- .check_validity_of_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
  
    if (is.null(features)) {
        features <- rownames(x)
    } else if (length(features) == 1) {
        features <- match.arg(features, c("type", "state", "none"))
        features <- get(paste0(features, "_markers"))(x)
        if (length(features) == 0)
            stop("No features matched the specified marker class.")
    } else {
        # replace problematic characters
        features <- gsub("-", "_", features)
        features <- gsub(":", ".", features)
        stopifnot(features %in% rownames(x))
    }

    # order clusters according to hierarchical 
    # clustering on median feature expressions 
    ms <- t(.agg(x[features, ], "cluster_id", "median"))
    d <- dist(ms, method="euclidean")
    o <- hclust(d, method="average")$order
    
    # construct data.frame & reference for plotting
    es <- assay(x[features, ], "exprs")
    cd <- data.frame(t(es), colData(x))
    df <- melt(cd, id.vars=names(colData(x)),
        variable.name="antigen", value.name="expression")
    df$ref <- "no"
    ref <- df
    ref$cluster_id <- "ref"
    ref$ref <- "yes"
    df <- rbind(df, ref)
    
    # compute cluster frequencies
    fq <- tabulate(x$cluster_id) / ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    
    # reorder clusters
    df$cluster_id <- factor(df$cluster_id, 
        levels=rev(c("ref", levels(x$cluster_id)[o])),
        labels=rev(c("ref", paste0(names(fq), " (", fq, "%)")[o])))
    
    ggplot(df, 
        aes_string(x="expression", y="cluster_id", col="ref", fill="ref")) + 
        facet_wrap(~antigen, scales="free_x", nrow=2) + 
        geom_density_ridges(alpha=.2) + 
        theme_ridges() + theme(
            legend.position="none",
            strip.background=element_blank(),
            strip.text=element_text(face="bold"))
}