# ==============================================================================
# Manual cluster merging
# ------------------------------------------------------------------------------

#' @rdname mergeClusters
#' @title Manual cluster merging
#'
#' @param x a \code{\link{daFrame}}.
#' @param table the merging table; a data.frame with columns 
#' \code{'old_cluster'}, \code{'new_cluster'} and \code{'label'}.
#' @param label a character string.
#' 
#' @return Writes the newly assigend cluster codes into the metadata slot 
#' \code{cluster_codes} of the input \code{daFrame} and returns the latter.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md, merging_table)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- mergeClusters(re, merging_table, "merging_1")
#' plotClusterHeatmap(re, k=20, m="merging_1", functional="CD4")
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

setMethod(f="mergeClusters", 
    signature=signature(x="daFrame"), 
    definition=function(x, table, label) {
        k <- max(table$old_cluster)
        m <- match(cluster_codes(x)[, k], table$old_cluster)
        new_ids <- table$new_cluster[m]
        metadata(x)$cluster_codes[, label] <- new_ids
        return(x)
    }
)
