#' @rdname extractClusters
#' @title Extract clusters from a \code{daFrame}
#' 
#' @description 
#' Extracts clusters from a \code{daFrame}. Populations will be 
#' either returned as a \code{flowSet} or written to FCS files,
#' depending on argument \code{as}.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
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
#' @author Mark D Robinson,
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @return a \code{flowSet} or character vector of the output file names.
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
#' extractClusters(re, k = "merging_1")
#' 
#' @importFrom flowCore write.FCS
# ------------------------------------------------------------------------------

setMethod(f="extractClusters", 
    signature=signature(x="daFrame"), 
    definition=function(x, k, clusters=NULL, 
        as=c("flowSet", "fcs"), out_dir=".", verbose=TRUE) {
        
        # check that cluster() has been run
        stopifnot("cluster_codes" %in% names(metadata(x)))
        stopifnot("cluster_id" %in% names(rowData(x)))
        
        # argument checks
        as <- match.arg(as)
        stopifnot(dir.exists(out_dir))
        stopifnot(is.logical(verbose))
        
        # get cluster IDs
        if (verbose) message("Extracting cluster IDs...")
        k <- check_validity_of_k(x, k)
        cluster_ids <- cluster_codes(x)[cluster_ids(x), k]
        if (is.null(clusters)) {
            clusters <- levels(cluster_ids)
        } else {
            stopifnot(all(clusters %in% levels(cluster_ids)))
        }
        
        # subset daFrame
        n <- sum(n_cells(x))
        x <- x[cluster_ids %in% clusters, ]
        if (verbose) message("Keeping ", nrow(x), "/", n, " cells...")
        
        # back-transform expressions
        cf <- metadata(x)$cofactor
        if (verbose) message("Back-transforming using cofactor ", cf, "...")
        es <- sinh(exprs(x))*cf
        
        inds <- split(seq_len(nrow(es)), cluster_ids)
        fct <- switch(as,
            flowSet = function(u, v) {
                ff <- new("flowFrame", exprs=es[u, , drop=FALSE])
                identifier(ff) <- v
                description(ff)$note <- paste0(nrow(ff), "/", n, " cells")
                if (verbose) message("Writing ", nrow(ff), " cells to flowFrame ", dQuote(v), "...")
                return(ff)
            },
            fcs = function(u, v) {
                ff <- new("flowFrame", exprs=es[u, , drop=FALSE])
                fn <- file.path(out_dir, paste0(v, ".fcs"))
                if (verbose) message("Writing ", nrow(ff), " cells to ", dQuote(fn), "...")
                suppressWarnings(write.FCS(ff, fn))
            }
        )
        ffs <- mapply(fct, inds, names(inds)) 
        if (as == "flowSet") as(ffs, "flowSet") else if (as == "fcs") ffs
    }
)
        