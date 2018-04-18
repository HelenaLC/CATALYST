# ==============================================================================
# FlowSOM clustering and ConsensusClusterPlus metaclustering
# ------------------------------------------------------------------------------
#' @rdname cluster
#' @title FlowSOM clustering and ConsensusClusterPlus metaclustering
#' 
#' @description 
#' \code{CATALYST::cluster()} runs \code{\link{FlowSOM}} clustering into 
#' \code{xdim}x\code{ydim} clusters, and \code{\link{ConsensusClusterPlus}} 
#' metaclustering into 2-\code{maxK} clusters. In the returned \code{daFrame}, 
#' those antigens used for clustering will be labelled as 'cell_type' markers, 
#' and the remainder of antigens as 'cell_state' markers.
#'
#' @param x a \code{\link{daFrame}}.
#' @param cols_to_use a character vector.
#' Specifies which antigens to use for clustering.
#' @param xdim,ydim numerical values specifying the grid size of the
#' self-orginizing map. The default 10x10 grid will yield 100 clusters. 
#' @param maxK numerical value. Specifies the maximum 
#' number of clusters to evaluate in the metaclustering.
#' 
#' @return 
#' The function will add information to the following slots 
#' of the input \code{daFrame} (and return it):
#' \itemize{
#' \item{\code{rowData}\describe{\itemize{
#' \item{\code{cluster_id}:
#' each cell's clustering ID as inferred by \code{FlowSOM}.}}
#' }}
#' \item{\code{colData}\describe{\itemize{
#' \item{\code{marker_class}: 
#' \code{"cell_type"} or \code{"cell_state"}. 
#' Specifies if an antigen has been used for clustering or not, respectively.}}
#' }}
#' \item{\code{metadata}\describe{\itemize{
#' \item{\code{SOM_codes}:
#' A table with dimensions K x (# cell type markers), 
#' where K = \code{xdim} x \code{ydim}. Contains the SOM codes.}
#' \item{\code{cluster_codes}:
#' A table with dimensions K x (\code{maxK} + 1). 
#' Contains the cluster codes for all metaclustering.}
#' \item{\code{delta_area}: 
#' A \code{\link{ggplot}} object. See above for details.}}
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
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # specify antigens to use for clustering
#' lineage <- c("CD3", "CD45", "CD4", "CD20", "CD33", 
#'     "CD123", "CD14", "IgM", "HLA_DR", "CD7")
#' (re <- cluster(re, cols_to_use=lineage))
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import ConsensusClusterPlus ggplot2
#' @importFrom flowCore flowFrame
#' @importFrom FlowSOM BuildSOM ReadInput
#' @importFrom graphics hist
#' @importFrom matrixStats colQuantiles
#' @importFrom reshape2 melt
#' @export
# ==============================================================================

setMethod(f="cluster",
    signature=signature(x="daFrame"),
    definition=function(x, cols_to_use, xdim=10, ydim=10, maxK=20) {
        
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
            maxK=maxK, reps=100, distance="euclidean", plot="pdf"))
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
        colData(x)$marker_class <- data.frame(
            row.names=colnames(exprs(x)),
            marker_class=factor(c("cell_state", "cell_type")[
                as.numeric(colnames(exprs(x)) %in% cols_to_use)+1],
                levels=levels(marker_classes(x))))
        metadata(x)$SOM_codes <- som$map$codes
        metadata(x)$cluster_codes <- cluster_codes
        metadata(x)$delta_area <- plot_delta_area(mc)
        return(x)
    }
)
