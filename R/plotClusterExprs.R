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
#' @importFrom dplyr %>% group_by_ summarise_at
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay colData rowData
#' @export

plotClusterExprs <- function(x, k="meta20", features=NULL) {
    
    stopifnot(is(x, "SingleCellExperiment"))
    k <- .check_validity_of_k(x, k)
    nms <- c("SOM_codes", "cluster_codes", "delta_area")
    stopifnot(!is.null(x$cluster_id), nms %in% names(metadata(x)))
   
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
    
    # calculate median features expressions
    es <- assay(x, "exprs")[features, , drop = FALSE]
    dat <- data.frame(t(es), colData(x))
    dat$cluster_id <- kids <- cluster_ids(x, k)
    meds <- dat %>% group_by_("cluster_id") %>% 
        summarise_at(features, funs(median))
    
    # hierarchical clustering
    d <- dist(meds[, features], method="euclidean")
    h <- hclust(d, method="average")
    o <- h$order
    
    # constrcut data.frame & reference for plotting
    dat <- melt(dat, id.vars=names(colData(x)),
        variable.name="antigen", value.name="expression")
    dat$ref <- "no"
    ref <- dat
    ref$cluster_id <- "ref"
    ref$ref <- "yes"
    df <- rbind(dat, ref)
    
    # compute cluster frequencies
    fq <- tabulate(kids) / ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(kids)
    
    # reorder clusters
    df$cluster_id <- factor(df$cluster_id, 
        levels=rev(c("ref", levels(kids)[o])),
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