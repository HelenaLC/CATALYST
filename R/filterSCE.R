#' @rdname filterSCE
#' @title `SingleCellExperiment` filtering
#' 
#' @description Filters cells/features from a \code{SingleCellExperiment} 
#'   using conditional statements a la \code{dplyr}.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param ... conditional statements separated by comma.
#'   Only rows/columns where the condition evaluates to TRUE are kept.
#' @param k numeric or character string. Specifies the clustering to extract 
#'   populations from. Must be one of \code{names(cluster_codes(x))}.
#'   Defaults to the 1st clustering available.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @return a \code{SingleCellExperiment}.
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md, merging_table)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' # one condition only, remove a single sample
#' filterSCE(sce, condition == "Ref", sample_id != "Ref1")
#' 
#' # keep only a subset of clusters
#' filterSCE(sce, cluster_id %in% c(7, 8, 18), k = "meta20")
#' 
#' @importFrom dplyr filter mutate_all select
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment 
#'   SummarizedExperiment assays colData rowData assays 
#' @export

filterSCE <- function(x, ..., k = NULL) {
    
    rd <- vapply(rowData(x), as.character, character(nrow(x)))
    cd <- vapply(colData(x), as.character, character(ncol(x)))
    rd <- data.frame(i=seq_len(nrow(x)), rd, 
        check.names=FALSE, stringsAsFactors=FALSE)
    cd <- data.frame(i=seq_len(ncol(x)), cd, 
        check.names=FALSE, stringsAsFactors=FALSE)
    
    # get cluster IDs for specified clustering
    if (is.null(k)) 
        k <- names(cluster_codes(x))[1] 
    k <- .check_validity_of_k(x, k)   
    cd$cluster_id <- cluster_ids(x, k)
    
    # filter rows & columns
    rdf <- try(dplyr::filter(rd, ...), silent=TRUE)
    cdf <- try(dplyr::filter(cd, ...), silent=TRUE)
    if (inherits(rdf, "try-error")) rdf <- rd
    if (inherits(cdf, "try-error")) cdf <- cd
    ri <- rdf$i; rdf <- select(rdf, -"i")
    ci <- cdf$i; cdf <- select(cdf, -"i")
    
    # convert to factors
    rdf <- mutate_all(rdf, factor)
    cdf <- mutate_all(cdf, factor)
    
    # update experimental design table
    md <- metadata(x)
    if (nrow(cdf) != ncol(x)) {
        ei <- metadata(x)$experiment_info
        cols <- intersect(colnames(cdf), colnames(ei))
        keep <- vapply(cols, function(u) 
            ei[, u] %in% levels(cdf[, u]), 
            logical(nrow(ei)))
        ei <- ei[apply(keep, 1, all), ]
        rownames(ei) <- NULL
        n_cells <- table(cdf$sample_id)
        m <- match(ei$sample_id, levels(cdf$sample_id))
        ei$n_cells <- as.numeric(n_cells[m])
        md$experiment_info <- ei
    }
    
    # revert colData(x)$cluster_id to 100 SOM clusters
    # and refactor colData factor columns
    cdf$cluster_id <- factor(x$cluster_id[ci], levels=levels(x$cluster_id))
    for (i in colnames(cdf)) if (i %in% names(ei))
        cdf[[i]] <- droplevels(factor(cdf[[i]], levels=levels(ei[[i]])))
    
    # subset reduced dimensions
    if (length(reducedDims(x)) > 0) {
        dr <- reducedDims(x)
        dr <- lapply(dr, "[", i=ci, j=TRUE)
    } else {
        dr <- NULL
    }
    
    # returned filtered SCE
    SingleCellExperiment(
        assays=lapply(assays(x), "[", i=ri, j=ci), 
        rowData=rdf, colData=cdf, 
        reducedDims=dr, metadata=md)
}