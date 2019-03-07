#' @rdname filter
#' @title Filter daFrame
#' 
#' @description 
#' Filters events/genes from a \code{daFrame} using conditional statements.
#'
#' @param x a \code{\link{daFrame}}.
#' @param ... conditional statements separated by comma.
#'   Only rows/columns where the condition evaluates to TRUE are kept.
#' @param k numeric or character string. Specifies the clustering to extract 
#'   populations from. Must be one of \code{names(cluster_codes(x))}.
#'   Defaults to \code{"som100"}.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
#' 
#' @return a \code{daFrame}.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md, merging_table)
#' daf <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run clustering
#' daf <- cluster(daf)
#' 
#' # one condition only, remove a single sample
#' filter(daf, condition == "Ref", sample_id != "Ref1")
#' 
#' # keep only a subset of clusters
#' filter(daf, cluster_id %in% c(7, 8, 18), k = "meta20")
#' 
#' @importFrom dplyr filter mutate_all select
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData rowData assays SummarizedExperiment
# ------------------------------------------------------------------------------

setMethod(f="filter", 
    signature=signature(x="daFrame", k = "character"), 
    definition=function(x, ..., k = NULL) {
        rd <- vapply(rowData(x), as.character, character(nrow(x)))
        cd <- vapply(colData(x), as.character, character(ncol(x)))
        rd <- data.frame(i=seq_len(nrow(x)), rd, 
            check.names=FALSE, stringsAsFactors=FALSE)
        cd <- data.frame(i=seq_len(ncol(x)), cd, 
            check.names=FALSE, stringsAsFactors=FALSE)

        # get cluster IDs for specified clustering
        k <- .check_validity_of_k(x, k)
        rd$cluster_id <- .get_cluster_ids(x, k)
        
        # filter rows & columns
        rdf <- try(filter(rd, ...), silent=TRUE)
        cdf <- try(filter(cd, ...), silent=TRUE)
        if (inherits(rdf, "try-error")) rdf <- rd
        if (inherits(cdf, "try-error")) cdf <- cd
        ri <- rdf$i; rdf <- select(rdf, -"i")
        ci <- cdf$i; cdf <- select(cdf, -"i")
        
        # convert to factors
        rdf <- mutate_all(rdf, factor)
        cdf <- mutate_all(cdf, factor)
        # fix marker_class levels
        cdf$marker_class <- factor(cdf$marker_class, 
            levels=levels(marker_classes(x))) 
        
        # update experimental design table
        ei <- metadata(x)$experiment_info
        cols <- intersect(colnames(rdf), colnames(ei))
        keep <- vapply(cols, function(u) 
            ei[, u] %in% levels(rdf[, u]), 
            logical(nrow(ei)))
        ei <- ei[apply(keep, 1, all), ]
        rownames(ei) <- NULL

        md <- metadata(x)
        if (nrow(x) != nrow(rdf)) {
            # update metadata
            md$experiment_info <- ei
            md$n_cells <- table(rdf$sample_id)
            md$tsne$Y <- md$tsne$Y[md$tsne_inds %in% ri, ]
            inds <- match(md$tsne_inds, ri)
            md$tsne_inds <- inds[!is.na(inds)]
        }
        
        # revert colData(x)$cluster_id to 100 SOM clusters
        rdf$cluster_id <- factor(cluster_ids(x)[ri], 
            levels=levels(cluster_ids(x)))
        
        # returned filtered daFrame
        se <- SummarizedExperiment(
            assays=lapply(assays(x), "[", i=ri, j=ci), 
            rowData=rdf, colData=cdf, metadata=md)
        as(se, "daFrame")
    }
)

#' @rdname filter
setMethod(f="filter", 
    signature=signature(x="daFrame", k = "missing"), 
    definition=function(x, ..., k = NULL) {
        filter(x, ..., k = names(cluster_codes(x))[1])
    }
)