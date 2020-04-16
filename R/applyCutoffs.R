#' @rdname applyCutoffs
#' @title Single-cell debarcoding (2)
#' @description Applies separation and mahalanobies distance cutoffs.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay character string specifying which assay data to use.
#'   Should be one of \code{assayNames(x)} and 
#'   correspond to expression-like not count data.
#' @param mhl_cutoff numeric mahalanobis distance threshold above which events 
#'   should be unassigned; ignored if \code{metadata(x)$mhl_cutoff} exists.
#' @param sep_cutoffs non-negative numeric of length one or of same length 
#'   as the number of rows in the \code{bc_key(x)}. Specifies the distance 
#'   separation cutoffs between positive and negative barcode populations 
#'   below which events should be unassigned. If \code{NULL} (default), 
#'   \code{applyCutoffs} will try to access \code{metadata(x)$sep_cutoffs}.
#'
#' @return the input \code{SingleCellExperiment} \code{x} is returned with 
#' updated \code{colData} columns \code{"bc_id"} and \code{"mhl_dist"}, 
#' and an additional \code{int_metadata} slot \code{"mhl_cutoff"} 
#' containing the applied mahalanobies distance cutoff.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
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
#' sce <- prepData(sample_ff)
#'     
#' # assign preliminary barcode IDs
#' # & estimate separation cutoffs
#' sce <- assignPrelim(sce, sample_key)
#' sce <- estCutoffs(sce)
#' 
#' # use estimated population-specific 
#' # vs. global separation cutoff(s)
#' sce1 <- applyCutoffs(sce)
#' sce2 <- applyCutoffs(sce, sep_cutoffs = 0.35)
#' 
#' # compare yields after applying cutoff(s)
#' c(global = mean(sce1$bc_id != 0), 
#' specific = mean(sce2$bc_id != 0))
#'   
#' @importFrom Matrix colMeans solve
#' @importFrom stats cov mahalanobis
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay assayNames
#' @export

applyCutoffs <- function(x, assay = "exprs", 
    mhl_cutoff = 30, sep_cutoffs = NULL) {
    # check validity of input arguments
    args <- as.list(environment())
    .check_args_applyCutoffs(args)
    
    chs <- channels(x)
    ms <- .get_ms_from_chs(chs)
    bc_key <- metadata(x)$bc_key
    bc_cols <- match(colnames(bc_key), ms)
    names(ids) <- ids <- rownames(bc_key)
    
    if (is.null(sep_cutoffs)) {
        seps <- metadata(x)$sep_cutoffs
        seps[is.na(seps)] <- 1
    } else {
        n_ids <- length(ids)
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
    }
    # split cells by barcode population
    cs <- split(seq_len(ncol(x)), x$bc_id)[ids]
    ns <- vapply(cs, length, numeric(1))
    # exclude events falling below separation cutoff
    ex_sep <- lapply(ids[ns != 0], function(id)
        x$delta[cs[[id]]] < seps[id])
    
    # compute mahalanobis distances given current separation cutoff
    n_bcs <- ncol(bc_key)
    mhl_dists <- rep(NA, ncol(x))
    for (id in ids[ns != 0]) {
        i <- cs[[id]][!ex_sep[[id]]]
        y <- assay(x, assay)[rowData(x)$is_bc, i, drop = FALSE]
        y <- t(as.matrix(y))
        if (length(y) != n_bcs && nrow(y) > n_bcs) {
            cvm <- cov(y)
            test <- tryCatch(
                solve(cvm) %*% cvm, 
                error = function(e) e)
            if (!inherits(test, "error")) 
                mhl_dists[i] <- mahalanobis(y, colMeans(y), cvm)
        }
    }
    ex_mhl <- mhl_dists > mhl_cutoff
    ex_mhl[is.na(ex_mhl)] <- FALSE
    x$bc_id[ex_mhl] <- x$bc_id[unlist(cs)[unlist(ex_sep)]] <- 0
    
    x$mhl_dist <- mhl_dists
    metadata(x)$mhl_cutoff <- mhl_cutoff
    return(x)
}