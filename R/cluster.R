# ==============================================================================
# tSNE-method for dbFrame-class
# ------------------------------------------------------------------------------
#' @rdname cluster
#' @title MDS plot
#' 
#' @description 
#' Multi-dimensional scaling (MDS) plot on median marker expressions.
#'
#' @param x a \code{\link{daFrame}}.
#' @param color numeric value between 2 and 20 OR a character string specifying
#' an antibody that appears in the metadata table of the input \code{daFrame}.
#' Specifies the color coding.
#' @param facette one of \code{NULL}, \code{"sample_id"} or \code{"condition"}.
#' Specifies whether the data should be split between sample IDs or conditions,
#' respectively. Defaults to NULL.
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- tSNE(re)
#' plotSNE(re, "CD4")
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import ggplot2
# ==============================================================================

setMethod(f="cluster",
    signature=signature(x="daFrame"),
    definition=function(x, cols_to_use, facette=NULL) {
        
        # flowSOM clustering
        message("o running FlowSOM clustering...")
        fsom <- ReadInput(flowFrame(exprs(x)))
        som <- BuildSOM(fsom, colsToUse=cols_to_use, silent=TRUE)
        codes <- som$map$codes
        cluster_ids <- som$map$mapping[, 1]
        
        # metaclustering
        message("o running ConsensusClusterPlus metaclustering...")
        pdf(NULL)
        mc <- suppressMessages(ConsensusClusterPlus(t(codes), 
            maxK=20, reps=100, distance="euclidean", plot="pdf"))
        dev.off()

        # get cluster codes
        cluster_codes <- data.frame(matrix(0, 100, 20, 
            dimnames=list(NULL, c(100, 2:20))), check.names=FALSE)
        for (k in seq_len(20)[-1])
            cluster_codes[, k] <- mc[[k]]$consensusClass
        
        col_data <- data.frame(row.names=colnames(exprs(x)),
            lineage=as.numeric(colnames(exprs(x)) %in% cols_to_use))
        col_data$functional <- as.numeric(!col_data$lineage)
        
        cluster_codes[, 1] <- seq_len(100)
        rowData(x)$cluster_id <- cluster_ids
        colData(x) <- col_data
        metadata(x)$cluster_codes <- cluster_codes
        
        
    }
)