# ==============================================================================
# Compute compensation matrix
# ------------------------------------------------------------------------------

#' @rdname computeSpillmat
#' @title Compute spillover matrix
#' 
#' @description 
#' Computes a spillover matrix from a priori 
#' identified single-positive populations.
#'
#' @param x 
#' a \code{\link{dbFrame}}.
#' @param method 
#' function to be used for computing spillover estimates
#' (see below for details).
#' @param interactions
#' \code{"default"} or \code{"all"}. Specifies which interactions spillover 
#' should be estimated for. The default exclusively takes into consideration 
#' interactions that are sensible from a chemical and physical point of view
#' (see below for more details).
#' @param trim
#' trim value used for estimation of spill values. 
#' Note that \code{trim = 0.5} is equivalent to using medians.
#' @param th
#' a single non-negative numeric. Specifies a threshold value below which spill
#' estimates will be set to 0. Applies only if \code{interaction="all"}.
#'
#' @return
#' Returns a square compensation matrix with dimensions and dimension names 
#' matching those of the input flowFrame. Spillover is assumed to be linear,
#' and, on the basis of their additive nature, spillover values are computed 
#' independently for each interacting pair of channels. 
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
#' M+/-1 (detection sensitivity), same metals (isotopic impurites) and M+16 
#' (oxide formation). By default, diagonal entries are set to 1. 
#' 
#' \code{interaction="all"} will estimate spill for all n x n - n 
#' interactions, where n denotes the number of single-color controls 
#' (= \code{nrow(bc_key(re))}).
#' 
#' @examples
#' # get single-stained control samples
#' data(ss_exp)
#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#' # debarcode single-positive populations
#' re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' head(computeSpillmat(x = re))
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom stats median

# ------------------------------------------------------------------------------

setMethod(f="computeSpillmat", 
    signature=signature(x="dbFrame"), 
    definition=function(x, method="default", interactions="default", 
        trim = .5, th = 10e-6) {
        
        if (sum(rowSums(bc_key(x)) == 1) != ncol(bc_key(x))) 
            stop("Cannot compute spillover matrix 
                from non single-staining experiment.")
        
        # check validity of input arguments
        if (!method %in% c("default", "classic"))
            stop("Invalid 'method' specified.\n", 
                "Valid options are \"default\" and \"classic\".\n",
                "See ?computeSpillmat for more details.")
        if (!interactions %in% c("default", "all"))
            stop("Invalid 'interactions' specified.\n", 
                "Valid options are \"default\" and \"all\".\n",
                "See ?computeSpillmat for more details.")
        
        # get intensities, no. of channels, masses and metals
        es <- exprs(x)
        chs <- colnames(es)
        ms <- as.numeric(regmatches(chs, gregexpr("[0-9]+", chs)))
        mets <- gsub("[[:digit:]]+Di", "", chs)
        
        # get barcode IDs and barcodes masses
        ids <- unique(bc_ids(x))
        ids <- ids[ids != 0]
        bc_ms <- as.numeric(rownames(bc_key(x)))
        
        # find which columns of loaded FCS file 
        # correspond to masses listed in barcode key
        bc_cols <- vapply(bc_ms, function(x) which(ms == x), numeric(1))
        
        if (interactions == "default") {
            # for each channel, get spillover candidate channels
            # (+/-1M, -16M and channels measuring isotopes)
            spill_cols <- get_spill_cols(ms, mets)
        } else if (interactions == "all") {
            # consider all channels
            spill_cols <- lapply(ms, function(x) which(ms != x & !is.na(ms)))
        }
        
        # compute and return compensation matrix
        SM <- diag(ncol(es))
        if (method == "default") {
            for (id in ids) {
                i <- bc_cols[bc_ms == id]
                j <- spill_cols[[which(ms == id)]]
                pos <- bc_ids(x) == id
                neg <- !bc_ids(x) %in% c(0, id, ms[spill_cols[[i]]])
                if (sum(neg) != 0) {
                    bg <- median(es[neg, j]) / median(es[neg, i])
                    if (is.na(bg)) 
                        bg <- 0
                } else {
                    bg <- 0
                }
                s <- es[pos, j] / es[pos, i] - bg
                s <- matrix(s, ncol=length(j))
                s[is.na(s) | s < 0] <- 0
                SM[i, j] <- matrixStats::colMedians(s)
            }
        } else {
            for (id in ids) {
                i <- bc_cols[bc_ms == id]
                j <- spill_cols[[which(ms == id)]]
                pos <- bc_ids(x) == id
                neg <- !bc_ids(x) %in% c(0, id, ms[spill_cols[[i]]])
                if (sum(neg) != 0) {
                    s <- (apply(es[pos, j], 2, mean, trim) - 
                            apply(es[neg, j], 2, mean, trim)) /
                        (mean(es[pos, i], trim) - mean(es[neg, i], trim))
                } else {
                    s <- apply(es[pos, j], 2, mean, trim) / 
                        mean(es[pos, i], trim)
                }
                s[is.na(s) | s < 0] <- 0
                SM[i, j] <- s
            }
        }
        colnames(SM) <- rownames(SM) <- chs
        if (interactions == "all") 
            SM[SM < th] <- 0
        SM[bc_cols, !is.na(ms)]
    })