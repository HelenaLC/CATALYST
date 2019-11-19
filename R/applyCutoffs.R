#' @rdname applyCutoffs
#' @title Single-cell debarcoding (2)
#' @description Applies separation and mahalanobies distance cutoffs.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param mhl_cutoff 
#'   mahalanobis distance threshold above which events should be unassigned. 
#'   This argument will be ignored if the \code{mhl_cutoff} slot of the input 
#'   \code{dbFrame} is specified.
#' @param sep_cutoffs 
#'   non-negative numeric of length one or of same length as the number of rows 
#'   in the \code{bc_key(x)}. Specifies the distance separation cutoffs between 
#'   positive and negative barcode populations below which events should be 
#'   unassigned. If \code{NULL} (default), \code{applyCutoffs} will try to 
#'   access the \code{sep_cutoffs} slot of the input \code{dbFrame}.
#'
#' @return the input \code{SingleCellExperiment} \code{x} is returned with 
#' updated \code{colData} columns \code{"bc_id"} and \code{"mhl_dist"}, 
#' and an additional \code{int_metadata} slot \code{"mhl_cutoff"} 
#' containing the applied mahalanobies distance cutoff.
#' 
#' @author Helena L. Crowell
#' 
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @examples
#' library(SingleCellExperiment)
#' 
#' # construct SCE
#' data(sample_ff, sample_key)
#' es <- as.matrix(exprs(sample_ff))
#' sce <- SingleCellExperiment(
#'     assays = list(counts = t(es)),
#'     rowData = pData(parameters(sample_ff)))
#'     
#' # assign preliminary barcode IDs
#' sce <- assignPrelim(x = sce, bc_key = sample_key)
#' 
#' # estimate population-specific cutoffs
#' sce2 <- estCutoffs(x = sce)
#' sce <- applyCutoffs(x = sce2, sep_cutoffs = 0.6)
#' plotEvents(sce, "A1", out_path = "~/Desktop", n = 1e3)
#'
#' @importFrom Matrix rowMeans solve
#' @importFrom methods is
#' @importFrom stats cov mahalanobis
#' @importFrom SingleCellExperiment int_metadata int_metadata<-
#' @importFrom SummarizedExperiment assay assayNames
#' @export

applyCutoffs <- function(x, 
    altExp = "barcodes", assay = "exprs", 
    mhl_cutoff = 30, sep_cutoffs = NULL) {
    stopifnot(
        is(x, "SingleCellExperiment"),
        is.character(assay), length(assay) == 1,
        is.null(altExp) || is.character(altExp) && 
            length(altExp) == 1 && altExp %in% altExpNames(x),
        is.numeric(mhl_cutoff), length(mhl_cutoff) == 1)
    
    if (is.null(altExp)) y <- x else y <- altExp(x, altExp)
    stopifnot(assay %in% assayNames(y),
        !is.null(y$bc_id), !is.null(y$delta),
        is.null(sep_cutoffs) | int_metadata(y)$sep_cutoffs)
    
    bc_key <- metadata(y)$bc_key
    ms <- .get_ms_from_chs(rownames(y))
    bc_cols <- match(colnames(bc_key), ms)
    n_ids <- length(names(ids) <- ids <- rownames(bc_key))
    
    if (!is.null(sep_cutoffs)) {
        seps <- sep_cutoffs
        if (length(seps) == 1) {
            seps <- rep(seps, n_ids)
            names(seps) <- ids
        } else {
            stopifnot(length(seps) == n_ids)
            if (is.null(names(seps))) {
                names(seps) <- ids
            } else {
                stopifnot(!all(names(seps) %in% ids))
            }
        }
    } else {
        seps <- int_metadata(y)$sep_cutoffs
        seps[is.na(seps)] <- 1
    }
    
    # compute mahalanobis distances given current separation cutoff
    cs <- split(seq_len(ncol(y)), y$bc_id)
    ns <- vapply(cs, length, numeric(1))
    ex_sep <- lapply(ids, function(id) y$delta[cs[[id]]] < seps[id])
    mhl_dists <- lapply(ids, function(id) {
        if (ns[id] == 0) next
        z <- assay(y, assay)[, cs[[id]][!ex_sep[[id]]]]
        if (length(z) != n_bcs & ncol(z) > n_bcs) {
            # compute covarianve matrix
            cov_mat <- cov(t(z))
            # check that it's invertible
            test <- tryCatch(
                solve(cov_mat) %*% cov_mat, 
                error = function(e) e)
            if (!inherits(test, "error")) {
                mahalanobis(t(z), rowMeans(z), cov_mat)
            }
        } 
    })
    ex_sep <- unlist(ex_sep)
    ex_mhl <- unlist(mhl_dists) > mhl_cutoff
    y$mhl_dist <- unlist(mhl_dists)
    y$bc_id[ex_sep] <- y$bc_id[!ex_sep][ex_mhl] <- 0
    int_metadata(y)$mhl_cutoff <- mhl_cutoff
    if (!is.null(altExp)) altExp(x, altExp) <- y else x <- y
    return(x)
}