#' @rdname mergeClusters
#' @title Manual cluster merging
#'
#' @description \code{mergeClusters} provides a simple wrapper 
#' to store a manual merging inside the input \code{SingleCellExperiment}.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k character string specifying the clustering to merge;
#'   valid values are \code{names(cluster_codes(x))}.
#' @param table merging table with 2 columns containing the cluster IDs to
#'   merge in the 1st, and the cluster IDs to newly assign in the 2nd column.
#' @param id character string used as a label for the merging.
#' @param overwrite logical specifying whether to force overwriting
#'   should a clustering with name \code{id} already exist.
#' 
#' @details in the following code snippets, 
#' \code{x} is a \code{SingleCellExperiment} object.
#' \itemize{
#' \item{merging codes are accesible through \code{cluster_codes(x)$id}}
#' \item{all functions that ask for specification of a clustering 
#'   (e.g. \code{\link{plotAbundances}}, \code{\link{plotMultiHeatmap}})
#'   take the merging ID as a valid input argument.}}
#' 
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}} 
#' with newly added cluster codes stored in \code{cluster_codes(.)$id}.
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
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md, merging_table)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' # merge clusters
#' sce <- mergeClusters(sce, 
#'   k = "meta20", 
#'   id = "merging",
#'   table = merging_table)
#' 
#' # tabulate manual merging
#' table(cluster_ids(sce, k = "merging"))
#' 
#' # visualize median type-marker expression
#' plotExprHeatmap(sce, 
#'   features = "type", 
#'   by = "cluster_id", 
#'   k = "merging",
#'   bars = TRUE)
#' 
#' @importFrom methods is
#' @importFrom SingleCellExperiment colData
#' @importFrom S4Vectors metadata
#' @export

mergeClusters <- function(x, k, table, id, overwrite = FALSE) {
    # check validity of input arguments
    .check_sce(x, TRUE)
    .check_k(x, k)
    table <- data.frame(table)
    stopifnot(
        is.character(id), length(id) == 1, 
        dim(table) != 0, ncol(table) == 2,
        nrow(table) == length(unique(table[, 1])),
        all(table[, 1] %in% levels(cluster_codes(x)[, k])),
        is.logical(overwrite), length(overwrite) == 1)
    
    if (!overwrite && id %in% names(cluster_codes(x)))
        stop("There already exists a clustering ", dQuote(id), ";\n",
            "specify a different 'id' or use 'overwrite = TRUE'.")
    
    m <- match(cluster_codes(x)[, k], table[, 1])
    new_ids <- table[m, 2]
    metadata(x)$cluster_codes[, id] <- factor(new_ids)
    return(x)
}