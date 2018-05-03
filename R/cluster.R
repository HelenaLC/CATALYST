#' @rdname cluster
#' @title FlowSOM clustering & ConsensusClusterPlus metaclustering
#' 
#' @description 
#' \code{cluster} will first group cells into \code{xdim}x\code{ydim} 
#' clusters using \pkg{FlowSOM}, and subsequently perform metaclustering 
#' with \pkg{ConsensusClusterPlus} into 2 through \code{maxK} clusters. 
#' In the returned \code{daFrame}, those antigens used for clustering will be 
#' labelled as '\code{type}' markers, and the remainder of antigens as 
#' '\code{state}' markers.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param cols_to_use 
#'   a character vector. Specifies which antigens to use for clustering.
#' @param xdim,ydim 
#'   numeric. Specify the grid size of the self-orginizing map. 
#'   The default 10x10 grid will yield 100 clusters. 
#' @param maxK 
#'   numeric. Specifies the maximum number of clusters to evaluate
#'   in the metaclustering. For \code{maxK = 20}, for example, 
#'   metaclustering will be performed for 2 through 20 clusters.
#' @param seed
#'   numeric. Sets random seed in \code{ConsensusClusterPlus()}.
#' 
#' @return 
#' The function will add information to the following slots 
#' of the input \code{daFrame} (and return it):
#' \itemize{
#' \item{\code{rowData}\itemize{
#' \item{\code{cluster_id}:
#'   each cell's cluster ID as inferred by \code{FlowSOM}. 
#'   One of 1, ..., \code{xdim}x\code{ydim}.}
#' }}
#' \item{\code{colData}\itemize{
#' \item{\code{marker_class}: 
#'   \code{"type"} or \code{"state"}. 
#'   Specifies whether an antigen has been used for clustering 
#'   or not, respectively.}
#' }}
#' \item{\code{metadata}\itemize{
#' \item{\code{SOM_codes}:
#'   a table with dimensions K x (# cell type markers), 
#'   where K = \code{xdim} x \code{ydim}. Contains the SOM codes.}
#' \item{\code{cluster_codes}:
#'   a table with dimensions K x (\code{maxK} + 1). 
#'   Contains the cluster codes for all metaclustering.}
#' \item{\code{delta_area}: 
#'   a \code{\link{ggplot}} object. See above for details.}
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
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # specify antigens to use for clustering
#' lineage <- c("CD3", "CD45", "CD4", "CD20", "CD33", 
#'     "CD123", "CD14", "IgM", "HLA_DR", "CD7")
#' (re <- cluster(re, cols_to_use=lineage))
#' 
#' @import ConsensusClusterPlus ggplot2
#' @importFrom flowCore flowFrame
#' @importFrom FlowSOM BuildSOM ReadInput
#' @importFrom graphics hist
#' @importFrom matrixStats colQuantiles
#' @importFrom reshape2 melt
# ------------------------------------------------------------------------------

setMethod(f="cluster",
    signature=signature(x="daFrame"),
    definition=function(x, cols_to_use, xdim=10, ydim=10, maxK=20, seed=1) {
        
        # replace dash with underscore
        cols_to_use <- gsub("-", "_", cols_to_use)
        
        # flowSOM clustering
        message("o running FlowSOM clustering...")
        fsom <- ReadInput(flowFrame(exprs(x)))
        som <- BuildSOM(fsom, colsToUse=cols_to_use, 
            silent=TRUE, xdim=xdim, ydim=ydim)
        
        # metaclustering
        message("o running ConsensusClusterPlus metaclustering...")
        pdf(NULL)
        mc <- suppressMessages(ConsensusClusterPlus(t(som$map$codes), 
            maxK=maxK, reps=100, distance="euclidean", seed=seed, plot="pdf"))
        dev.off()
        
        # get cluster codes
        k <- xdim * ydim
        mcs <- seq_len(maxK)[-1]
        cluster_codes <- data.frame(factor(seq_len(k)), 
            sapply(mc[-1], function(x) factor(x$consensusClass)))
        colnames(cluster_codes) <- c(k, mcs)
        # reorder factor levels
        cluster_codes <- lapply(cluster_codes, function(codes) 
            factor(codes, levels=sort(as.numeric(levels(codes)))))
        cluster_codes <- data.frame(cluster_codes, check.names=FALSE)

        rowData(x)$cluster_id <- as.factor(som$map$mapping[, 1])
        colData(x)$marker_class <- factor(c("state", "type")[
            as.numeric(colnames(exprs(x)) %in% cols_to_use)+1],
            levels=levels(marker_classes(x)))
        metadata(x)$SOM_codes <- som$map$codes
        metadata(x)$cluster_codes <- cluster_codes
        metadata(x)$delta_area <- plot_delta_area(mc)
        return(x)
    }
)
