#' @rdname computeSpillmat
#' @title Compute spillover matrix
#' 
#' @description Computes a spillover matrix from 
#' priorly identified single-positive populations.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay character string specifying which assay to use; should be one
#'   of \code{assayNames(x)} and correspond to count-like data, as linearity 
#'   assumptions underlying spillover estimation won't hold otherwise.
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
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#'
#' @references
#' Coursey, J.S., Schwab, D.J., Tsai, J.J., Dragoset, R.A. (2015).
#' Atomic weights and isotopic compositions,
#' (available at http://physics.nist.gov/Comp).
#'
#' @examples
#' # construct SCE from single-stained control samples
#' data(ss_exp)
#' sce <- prepData(ss_exp)
#' 
#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#'
#' # debarcode single-positive populations
#' sce <- assignPrelim(sce, bc_ms)
#' sce <- estCutoffs(sce)
#' sce <- applyCutoffs(sce)
#' 
#' # estimate & extract spillover matrix 
#' sce <- computeSpillmat(sce)
#' 
#' library(SingleCellExperiment)
#' head(metadata(sce)$spillover_matrix)
#'
#' @importFrom methods is
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay assayNames
#' @export

computeSpillmat <- function(x, assay = "counts", 
    interactions = c("default", "all"), 
    method = c("default", "classic"),
    trim = 0.5, th = 10e-6) {

    interactions <- match.arg(interactions)
    method <- match.arg(method)

    stopifnot(
        is(x, "SingleCellExperiment"), .check_assay(x, assay),
        !is.null(metadata(x)$bc_key), !is.null(x$bc_id),
        is.numeric(trim), length(trim) == 1, !trim < 0, !trim > 0.5,
        is.numeric(th), length(th) == 1)

    n_bcs <- nrow(bc_key <- metadata(x)$bc_key)
    if (sum(rowSums(bc_key) == 1) != n_bcs)
        stop("Cannot compute spillover matrix",
            " from non single-staining experiment.")

    # get channel masses & metals
    chs <- channels(x)
    ms <- .get_ms_from_chs(chs)
    mets <- .get_mets_from_chs(chs)
    
    # get barcode IDs & barcode channels
    ids <- setdiff(unique(x$bc_id), 0)
    bc_chs <- match(bc_ms <- colnames(bc_key), ms)

    # for each channel, get spillover candidate channels
    # (by default, +/-1M, -16M and channels measuring isotopes)
    if (interactions == "default") {
        ex <- spill_chs <- .get_spill_chs(ms, mets)
    } else if (interactions == "all") {
        is <- seq_along(ms)
        spill_chs <- lapply(is, function(i) is[-i])
        ex <- .get_spill_chs(ms, mets)
    }
    
    # initialize spillover matrix
    sm <- diag(nrow(x))
    dimnames(sm) <- list(chs, chs)
    
    # split cells by barcode population
    cs <- split(seq_len(ncol(x)), x$bc_id)
    y <- assay(x, assay)
    for (id in ids) {
        i <- match(id, ms)
        pos <- cs[[id]]
        # neg. pop. = unassigned, pop. itself,
        # and events assigned to spill-affected pops.
        neg <- !x$bc_id %in% c(0, id, ms[ex[[i]]])
        pos_i <- y[i, pos]
        neg_i <- y[i, neg]
        for (j in spill_chs[[i]]) {
            pos_j <- y[j, pos]
            # further exclude events assigned to population
            # for which interaction is calculated,
            # or that are spill-affected by it
            neg_j <- y[j, neg & !x$bc_id %in% c(ms[j], ms[ex[[j]]])]
            sm[i, j] <- .get_sij(pos_i, neg_i, pos_j, neg_j, method, trim)
        }
    }
    sm[sm < th] <- 0
    sm <- sm[bc_chs, !is.na(ms)]
    metadata(x)$spillover_matrix <- sm
    return(x)
}
