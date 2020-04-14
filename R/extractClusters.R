#' @rdname extractClusters
#' @title Extract clusters from a \code{SingleCellExperiment}
#' 
#' @description 
#' Extracts clusters from a \code{SingleCellExperiment}. 
#' Populations will be either returned as a \code{flowSet} 
#' or written to FCS files, depending on argument \code{as}.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k 
#'   numeric or character string. 
#'   Specifies the clustering to extract populations from.
#'   Must be one of \code{names(cluster_codes(x))}.
#' @param clusters
#'   a character vector.
#'   Specifies which clusters to extract.
#'   \code{NULL} = all clusters.
#' @param as
#'   \code{"flowSet"} or \code{"fcs"}. 
#'   Specifies whether clusters should be return 
#'   as a \code{flowSet} or written to FCS files.
#' @param out_dir
#'   a character string.
#'   Specifies where FCS files should be writen to.
#'   Defaults to the working directory.
#' @param verbose
#'   logical.
#'   Should information on progress be reported?
#' 
#' @author Mark D Robinson & Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @return a \code{flowSet} or character vector of the output file names.
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md, merging_table)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' # merge clusters
#' sce <- mergeClusters(sce, k="meta20", table=merging_table, id="merging_1")
#' extractClusters(sce, k="merging_1", clusters=c("NK cells", "surface-"))
#' 
#' @importFrom flowCore identifier<- keyword<- write.FCS
#' @importFrom methods is
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay colData rowData
#' @export

extractClusters <- function(x, k, clusters = NULL, 
    as = c("flowSet", "fcs"), out_dir = ".", verbose = TRUE) {
    
    # check validity of input arguments
    as <- match.arg(as)
    .check_sce(x, TRUE)
    k <- .check_k(x, k)
    stopifnot(
        is.character(out_dir), dir.exists(out_dir), 
        is.logical(verbose), length(verbose) == 1)
    
    # get cluster IDs
    if (verbose) message("Extracting cluster IDs...")
    cluster_ids <- cluster_ids(x, k)
    if (is.null(clusters)) {
        clusters <- levels(cluster_ids)
    } else {
        stopifnot(clusters %in% levels(cluster_ids))
    }
    
    # subset SCE
    n <- ncol(x)
    x <- x[, cluster_ids %in% clusters]
    if (verbose) message("Keeping ", ncol(x), "/", n, " cells...")
    
    # back-transform expressions
    es <- assay(x, "exprs")
    if (!is.null(cf <- int_metadata(x)$cofactor)) {
        if (verbose) message("Back-transforming using cofactor ", cf, "...")
        es <- sinh(es)*cf
    }
    
    fun <- switch(as,
        flowSet = function(u, v) {
            ff <- new("flowFrame", exprs=t(es[, u, drop=FALSE]))
            identifier(ff) <- v
            keyword(ff)$note <- paste0(nrow(ff), "/", n, " cells")
            if (verbose) message("Storing ", nrow(ff), 
                " cells in flowFrame ", dQuote(v), "...")
            return(ff)
        },
        fcs = function(u, v) {
            ff <- new("flowFrame", exprs=t(es[, u, drop=FALSE]))
            fn <- file.path(out_dir, paste0(v, ".fcs"))
            if (verbose) message("Writing ", nrow(ff), 
                " cells to ", dQuote(fn), "...")
            suppressWarnings(write.FCS(ff, fn))
        }
    )
    cs_by_k <- split(seq_len(ncol(es)), cluster_ids(x, k))
    ffs <- mapply(fun, cs_by_k, names(cs_by_k))[clusters]
    if (as == "flowSet") as(ffs, "flowSet") else if (as == "fcs") ffs
}