#' @rdname cluster
#' @title \code{FlowSOM} clustering & 
#' \code{ConsensusClusterPlus} metaclustering
#' 
#' @description 
#' \code{cluster} will first group cells into \code{xdim}x\code{ydim} 
#' clusters using \pkg{FlowSOM}, and subsequently perform metaclustering 
#' with \pkg{ConsensusClusterPlus} into 2 through \code{maxK} clusters.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param features a character vector specifying 
#'   which features to use for clustering; valid values are
#'   \code{"type"/"state"} for \code{type/state_markers(x)} 
#'   if \code{rowData(x)$marker_class} have been specified; 
#'   a subset of \code{rownames(x)}; NULL to use all features.
#' @param xdim,ydim numeric specifying the grid size of the 
#'   self-orginizing map; passed to \code{\link[FlowSOM]{BuildSOM}}. 
#'   The default 10x10 grid will yield 100 clusters.
#' @param maxK numeric specifying the maximum number of 
#'   clusters to evaluate in the metaclustering; passed to 
#'   \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}.
#'   The default (\code{maxK = 20}) will yield 2 through 20 metaclusters.
#' @param verbose logical. Should information on progress be reported?
#' @param seed numeric. Sets the random seed for reproducible results 
#'   in \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}.
#' 
#' @return a \code{SingleCellEcperiment} with the following newly added data:
#' \itemize{
#' \item{\code{colData}\itemize{
#' \item{\code{cluster_id}:
#'   each cell's cluster ID as inferred by \code{FlowSOM}. 
#'   One of 1, ..., \code{xdim}x\code{ydim}.}
#' }}
#' \item{\code{rowData}\itemize{
#' \item{\code{marker_class}: added when previosly unspecified. \code{"type"}
#'   when an antigen has been used for clustering, otherwise \code{"state"}.}
#' \item{\code{used_for_clustering}: logical indicating 
#'   whether an antigen has been used for clustering.}
#' }}
#' \item{\code{metadata}\itemize{
#' \item{\code{SOM_codes}:
#'   a table with dimensions K x (# cell type markers), 
#'   where K = \code{xdim} x \code{ydim}. Contains the SOM codes.}
#' \item{\code{cluster_codes}:
#'   a table with dimensions K x (\code{maxK} + 1). 
#'   Contains the cluster codes for all metaclustering.}
#' \item{\code{delta_area}: 
#'   a \code{\link{ggplot}} object (see details).}
#' }}
#' }
#' 
#' @details 
#' The delta area represents the amount of extra cluster stability gained when 
#' clustering into k groups as compared to k-1 groups. It can be expected that 
#' high stability of clusters can be reached when clustering into the number of 
#' groups that best fits the data. The "natural" number of clusters present in 
#' the data should thus corresponds to the value of k where there is no longer 
#' a considerable increase in stability (pleateau onset).
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
#' # construct SCE
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run clustering
#' (sce <- cluster(sce))
#' 
#' # view all available clustering
#' names(cluster_codes(sce))
#' 
#' # access specific clustering resolution
#' table(cluster_ids(sce, "meta8"))
#' 
#' # view delta area plot
#' delta_area(sce)
#' 
#' @import ConsensusClusterPlus ggplot2
#' @importFrom dplyr %>% mutate_all
#' @importFrom flowCore flowFrame
#' @importFrom FlowSOM BuildSOM ReadInput
#' @importFrom graphics hist
#' @importFrom matrixStats colQuantiles
#' @importFrom methods is
#' @importFrom purrr map
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment assay rowData rowData<-
#' @importFrom S4Vectors DataFrame metadata<-
#' @export

cluster <- function(x, features = "type",
    xdim = 10, ydim = 10, maxK = 20, verbose = TRUE, seed = 1) {
    
    # validity checks
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(is.logical(verbose), length(verbose) == 1,
        vapply(list(xdim, ydim, maxK, seed), function(arg) 
            is.numeric(arg) && length(arg) == 1, logical(1)))
    features <- .get_features(x, features)
    
    # assign marker classes if unspecified
    if (is.null(marker_classes(x))) {
        rowData(x)$marker_class <- factor(c("state", "type")[
            as.numeric(rownames(x) %in% features)+1], 
            levels = c("type", "state", "none"))
    }
    # store which markers were used for clustering
    rowData(x)$used_for_clustering <- rownames(x) %in% features
    
    # do 'FlowSOM' clustering
    if (verbose)
        message("o running FlowSOM clustering...")
    fsom <- ReadInput(flowFrame(t(assay(x, "exprs"))))
    som <- BuildSOM(fsom, colsToUse=features, 
        silent=TRUE, xdim=xdim, ydim=ydim)
    
    # do 'ConsensusClusterPlus' metaclustering
    if (verbose)
        message("o running ConsensusClusterPlus metaclustering...")
    pdf(NULL)
    mc <- suppressWarnings(suppressMessages(
        ConsensusClusterPlus(t(som$map$codes), 
            maxK = maxK, reps = 100, 
            distance = "euclidean", 
            seed = seed, plot = NULL)))
    dev.off()
    
    # get cluster codes
    k <- xdim * ydim
    mcs <- seq_len(maxK)[-1]
    
    # construct data.frame of clustering codes
    codes <- data.frame(seq_len(k), map(mc[-1], "consensusClass"))
    codes <- mutate_all(codes, function(u) factor(u, levels = sort(unique(u))))
    colnames(codes) <- c(sprintf("som%s", k), sprintf("meta%s", mcs))
    
    # store cluster assignments in cell metadata;
    # cluster / SOM codes & delta area plot in metadata
    x$cluster_id <- factor(som$map$mapping[, 1])
    metadata(x)$cluster_codes <- codes
    metadata(x)$SOM_codes <- som$map$codes
    metadata(x)$delta_area <- .plot_delta_area(mc)
    return(x)
}
