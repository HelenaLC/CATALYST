#' @rdname computeSpillmat
#' @title Compute spillover matrix
#' 
#' @description Computes a spillover matrix from 
#' priorly identified single-positive populations.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param method
#'   \code{"default"} or \code{"classic"}. Specifies the function
#'   to be used for spillover estimation (see below for details).
#' @param interactions
#'   \code{"default"} or \code{"all"}. Specifies which interactions spillover
#'   should be estimated for. The default exclusively takes into consideration
#'   interactions that are sensible from a chemical and physical point of view
#'   (see below for more details).
#' @param trim
#'   numeric. Specifies the trim value used for estimation of spill values.
#'   Note that \code{trim = 0.5} is equivalent to using medians.
#' @param th
#'   single non-negative numeric. Specifies the threshold value
#'   below which spill estimates will be set to 0.
#'
#' @details
#' The \code{default} method estimates the spillover as the median ratio
#' between the unstained spillover receiving and the stained spillover
#' emitting channel in the corresponding single stained populations.
#'
#' \code{method = "classic"} will compute the slope of a line through
#' the medians (or trimmed means) of stained and unstained populations.
#' The medians (or trimmed means) computed from events that are i) negative
#' in the respective channels; and, ii) not assigned to interacting channels;
#' and, iii) not unassigned are subtracted as to account for background.
#'
#' \code{interactions="default"} considers only expected interactions, that is,
#' M+/-1 (detection sensitivity), M+16 (oxide formation) and channels measuring
#' metals that are potentially contaminated by isotopic impurites
#' (see reference below and \code{\link{isotope_list}}).
#'
#' \code{interaction="all"} will estimate spill for all n x n - n
#' interactions, where n denotes the number of single-color controls
#' (= \code{nrow(bc_key(re))}).
#'
#' @return
#' Returns a square compensation matrix with dimensions and dimension names
#' matching those of the input flowFrame. Spillover is assumed to be linear,
#' and, on the basis of their additive nature, spillover values are computed
#' independently for each interacting pair of channels.
#'
#' @author Helena L. Crowell
#'
#' @references
#' Coursey, J.S., Schwab, D.J., Tsai, J.J., Dragoset, R.A. (2015).
#' Atomic weights and isotopic compositions,
#' (available at http://physics.nist.gov/Comp).
#'
#' @examples
#' # get single-stained control samples
#' data(ss_exp)
#' es <- as.matrix(exprs(ss_exp))
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(
#'     assays = list(counts = t(es)),
#'     rowData = pData(parameters(ss_exp)))

#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#'
#' # debarcode single-positive populations
#' sce <- assignPrelim(x = sce, bc_key = bc_ms)
#' sce <- estCutoffs(x = sce)
#' sce <- applyCutoffs(x = sce)
#' head(computeSpillmat(x = re))
#'
#' @importFrom methods is
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay assayNames altExp altExpNames
#' @export

computeSpillmat <- function(x, altExp = "barcodes", assay = "exprs", 
    interactions = c("default", "all"), method = c("default", "classic"),
    trim = 0.5, th = 10e-6) {

    interactions <- match.arg(interactions)
    method <- match.arg(methods)

    stopifnot(
        is(x, "SingleCellExperiment"),
        is.null(altExp) || is.character(altExp) && 
            length(altExp) == 1 && altExp %in% altExpNames(x),
        is.character(assay), length(assay) == 1, 
        is.null(altExp) && assay %in% assayNames(x) ||
            !is.null(altExp) && assay %in% assayNames(altExp(x, altExp)),
        is.numeric(trim), length(trim) == 1, !trim < 0, !trim > 0.5,
        is.numeric(th), length(th) == 1)

    if (is.null(altExp)) y <- x else y <- altExp(x, altExp)
    n_bcs <- nrow(bc_key <- metadata(y)$bc_key)
    if (sum(rowSums(bc_key) == 1) != n_bcs)
        stop("Cannot compute spillover matrix",
            " from non single-staining experiment.")

    # get channel masses & metals
    ms <- .get_ms_from_chs(rownames(x))
    mets <- .get_mets_from_chs(rownames(x))
    
    # get barcode IDs & barcodes masses
    ids <- setdiff(unique(y$bc_id), 0)
    bc_cols <- match(bc_ms <- colnames(bc_key), ms)

    # for each channel, get spillover candidate channels
    # (by default, +/-1M, -16M and channels measuring isotopes)
    if (interactions == "default") {
        spill_cols <- .get_spill_cols(ms, mets)
        ex <- spill_cols
    } else if (interactions == "all") {
        spill_cols <- lapply(ms, function(u) setdiff(ms, c(u, NA)))
        ex <- .get_spill_cols(ms, mets)
    }
    
    # initialize spillover matrix
    sm <- diag(nrow(x))
    chs <- rownames(x)
    dimnames(sm) <- list(chs, chs)
    
    # split cells by barcode population
    cs <- split(seq_len(ncol(y)), y$bc_id)
    z <- assay(x, assay)
    for (id in ids) {
        i <- match(id, ms)
        pos <- cs[[id]]
        neg <- !y$bc_id %in% c(0, id, ms[ex[[i]]])
        pos_i <- z[i, pos]
        neg_i <- z[i, neg]
        for (j in spill_cols[[i]]) {
            pos_j <- z[j, pos]
            # further exclude events assigned to population
            # for which interaction is calculated
            neg_j <- z[j, y$bc_id[neg] != ms[j] 
                & !y$bc_id[neg] %in% ms[ex[[j]]]]
            sm[i, j] <- .get_sij(pos_i, neg_i, pos_j, neg_j, method, trim)
        }
    }
    sm[sm < th] <- 0
    sm[bc_cols, !is.na(ms)]
}