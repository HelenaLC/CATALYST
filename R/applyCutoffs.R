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
#' assay(sce, "exprs") <- asinh(assay(sce, "counts") / 10)
#'     
#' # assign preliminary barcode IDs
#' # & estimate separation cutoffs
#' sce <- assignPrelim(x = sce, bc_key = sample_key)
#' sce <- estCutoffs(sce)
#' 
#' # apply cutoffs
#' sce <- estCutoffs(x = sce)
#' sce <- applyCutoffs(x = sce, sep_cutoffs = 0)
#' plotEvents(sce, "A1", out_path = "~/Desktop", n = 1e3)
#'
#' @importFrom Matrix rowMeans solve
#' @importFrom methods is
#' @importFrom stats cov mahalanobis
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay assayNames
#' @export

applyCutoffs <- function(x, assay = "exprs", 
    mhl_cutoff = 30, sep_cutoffs = NULL) {
    stopifnot(
        is(x, "SingleCellExperiment"),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        !is.null(x$bc_id), !is.null(x$delta),
        !is.null(sep_cutoffs) | !is.null(metadata(x)$sep_cutoffs))
    
    n_bcs <- ncol(bc_key <- metadata(x)$bc_key)
    ms <- .get_ms_from_chs(rownames(x))
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
                stopifnot(all(names(seps) %in% ids))
            }
        }
    } else {
        seps <- metadata(x)$sep_cutoffs
        seps[is.na(seps)] <- 1
    }
    
    # compute mahalanobis distances given current separation cutoff
    cs <- split(seq_len(ncol(x)), x$bc_id)
    ns <- vapply(cs, length, numeric(1))
    ex_sep <- ex_mhl <- logical(ncol(x))
    ex_sep[unlist(cs[ids])] <- unlist(lapply(ids, 
        function(id) x$delta[cs[[id]]] < seps[id]))
    mhl_dists <- lapply(ids, function(id) {
        if (ns[id] == 0) return(NULL)
        i <- cs[[id]][!ex_sep[cs[[id]]]]
        foo <- rep(NA, length(i))
        if (length(i) == 0) return(foo)
        y <- assay(x, assay)[, i]
        if (length(y) != n_bcs & ncol(y) > n_bcs) {
            # compute covarianve matrix
            cov_mat <- cov(t(y))
            # check that it's invertible
            test <- tryCatch(
                solve(cov_mat) %*% cov_mat, 
                error = function(e) e)
            if (!inherits(test, "error")) { 
                return(mahalanobis(t(y), rowMeans(y), cov_mat))
            } else {
                foo
            }
        } else {
            foo
        } 
    })
    x$mhl_dist <- NA
    x$mhl_dist[which(!ex_sep[unlist(cs[ids])])] <- unlist(mhl_dists)
    ex_mhl <- x$mhl_dist > mhl_cutoff
    ex_mhl[is.na(ex_mhl)] <- FALSE
    x$bc_id[ex_sep | ex_mhl] <- 0
    metadata(x)$mhl_cutoff <- mhl_cutoff
    return(x)
}