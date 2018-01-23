# ==============================================================================
# Manual cluster merging
# ------------------------------------------------------------------------------

#' @rdname mergeClusters
#' @title Manual cluster merging
#' 
#' @description 
#'
#' @param x a \code{\link{daFrame}}.
#' @param table the merging table; a data.frame with columns 
#' \code{'old_cluster'}, \code{'new_cluster'} and \code{'label'}.
#' @param label a character string.
#' 
#' @return Writes the newly assigend cluster codes into the metadata slot 
#' \code{cluster_codes} of the input \code{daFrame} and returns the latter.
#'
#' @references 
#' Nowicka M, Krieg C, Weber LM et al.
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md, merging_table)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- mergeClusters(re, merging_table, "merging_1")
#' plotClusterHeatmap(re, k=20, m="merging_1", functional="CD4")
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2
# ==============================================================================

setMethod(f="mergeClusters", 
    signature=signature(x="daFrame"), 
    definition=function(x, table, label) {

        # check validity of input table
        table <- data.frame(table[order(table$old_cluster), ], row.names=NULL)
        k <- max(table$old_cluster)
        # if (!all(table$old_cluster %in% seq_len(k)))
        #     stop("Insufficient 'new_cluster' assignments provided.\n",
        #         "Please make sure 'old_cluster's 1-", k,
        #         " appear in the input 'table'.")
        unique_labels <- sapply(unique(table$new_cluster), function(id)
            length(unique(table[table$new_cluster == id, ]$label)) == 1)
        if (any(!unique_labels)) 
            stop("Invalid label(s) detected.\n       All labels provided", 
                " for the same 'new_cluster' must be identical.")
        
        m <- match(cluster_codes(x)[, k], table$old_cluster)
        new_ids <- table$label[m]
        metadata(x)$cluster_codes[, label] <- new_ids
        return(x)
    }
)