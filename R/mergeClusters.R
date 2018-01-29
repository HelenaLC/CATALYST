# ==============================================================================
# Manual cluster merging
# ------------------------------------------------------------------------------

#' @rdname mergeClusters
#' @title Manual cluster merging
#'
#' @description 
#' \code{mergeClusters} provides a simple wrapper to store a manual merging 
#' inside the input \code{daFrame}.
#'
#' @param x a \code{\link{daFrame}}.
#' @param table the merging table; a data.frame with columns 
#' \code{'old_cluster'}, \code{'new_cluster'} and \code{'label'}.
#' @param id a character string to use as a label for the merging.
#' 
#' @return Writes the newly assigend cluster codes into the metadata slot 
#' \code{cluster_codes} of the input \code{daFrame} and returns the latter.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md, merging_table)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run clustering
#' lineage <- c("CD3", "CD45", "CD4", "CD20", "CD33", 
#'     "CD123", "CD14", "IgM", "HLA_DR", "CD7")
#' re <- cluster(re, cols_to_use=lineage)
#' 
#' # merge clusters
#' re <- mergeClusters(re, merging_table, "merging_1")
#' plotClusterHeatmap(re, k=20, m="merging_1", type2="pNFkB")
#' 
#' @details 
#' in the following code snippets, \code{x} is a \code{daFrame} object.
#' \itemize{
#' \item merging codes are accesible through \code{cluster_codes(x)$id}
#' \item all functions that ask for specification of a clustering 
#' (e.g. \code{\link{plotAbundances}}, \code{\link{plotClusterHeatmap}})
#' take the merging ID as a valid input argument. 
#' }
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
    definition=function(x, table, id) {
        if (id %in% colnames(metadata(x)$cluster_codes)) {
            stop("There already exists a clustering named ",
                id, ". Please specify a different identifier.")
        }
        k <- max(table$old_cluster)
        m <- match(cluster_codes(x)[, k], table$old_cluster)
        new_ids <- table$new_cluster[m]
        metadata(x)$cluster_codes[, id] <- new_ids
        return(x)
    }
)
