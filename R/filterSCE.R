#' @rdname filterSCE
#' @title \code{SingleCellExperiment} filtering
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
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
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
    
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    
    # construct data.frames of cell & feature metadata
    rd <- data.frame(rowData(x), check.names = FALSE)
    cd <- data.frame(colData(x), check.names = FALSE)
    rd$i <- seq_len(nrow(x))
    cd$i <- seq_len(ncol(x))
    
    # get specified clustering IDs
    if (!is.null(k)) {
        stopifnot(
            !is.null(x$cluster_id),
            !is.null(cluster_codes(x)))
        k <- .check_k(x, k)   
        cd$cluster_id <- cluster_ids(x, k)
    }
    
    # filter rows & columns
    rdf <- try(dplyr::filter(rd, ...), silent = TRUE)
    cdf <- try(dplyr::filter(cd, ...), silent = TRUE)
    if (inherits(rdf, "try-error")) rdf <- rd
    if (inherits(cdf, "try-error")) cdf <- cd
    ri <- rdf$i; rdf <- select(rdf, -"i")
    ci <- cdf$i; cdf <- select(cdf, -"i")
    rdf <- droplevels(rdf)
    cdf <- droplevels(cdf)
    
    # update experimental design table
    md <- metadata(x)
    ei <- .get_ei(cdf)
    md$experiment_info <- ei
    
    # revert colData(x)$cluster_id to 100 SOM clusters
    if (!is.null(cluster_codes(x)))
        cdf$cluster_id <- factor(
            x$cluster_id[ci], 
            levels = levels(x$cluster_id))
    
    # refactor 'colData' factor columns
    for (i in colnames(cdf)) if (i %in% names(ei))
        cdf[[i]] <- droplevels(factor(cdf[[i]], levels=levels(ei[[i]])))
    
    # subset reduced dimensions
    if (length(reducedDims(x)) > 0) {
        dr <- reducedDims(x)
        dr <- lapply(dr, "[", i = ci, j = TRUE)
    } else dr <- NULL
    
    # subset assays & returned filtered SCE
    as <- lapply(assays(x), function(a) 
        a[ri, ci, drop = FALSE])
    SingleCellExperiment(assays = as, 
        rowData = rdf, colData = cdf, 
        reducedDims = dr, metadata = md)
}