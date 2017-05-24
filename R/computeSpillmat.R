# ==============================================================================
# Compute compensation matrix
# ------------------------------------------------------------------------------

#' @rdname computeSpillmat
#' @title Compute spillover matrx
#' 
#' @description 
#' Computes the spillover matrix based on 
#' a priori identified single-positive popultions.
#'
#' @param x 
#' a \code{\link{dbFrame}}.
#' @param method 
#' function to be used for computing spillover estimates. 
#' (see below for details)
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
#' The default method estimates the spillover as the median ratio between 
#' the unstained spillover receiving and the stained spillover emitting 
#' channel in the corresponding single stained populations. 
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
#' 
#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#' 
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
        
        # get no. of channels, masses and metals
        chs <- colnames(exprs(x))
        n_chs <- length(chs)
        ms <- as.numeric(regmatches(chs, gregexpr("[0-9]+", chs)))
        mets <- gsub("[[:digit:]]+Di", "", chs)
        
        # get barcode IDs and barcodes masses
        ids <- unique(bc_ids(x))
        ids <- sort(ids[which(ids != 0)])
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
            spill_cols <- lapply(ms, function(i) which(ms != i & !is.na(ms)))
        }
        
        # compute and return compensation matrix
        SM <- diag(n_chs)
        for (i in ids) {
            j <- bc_cols[ids == i]
            pos <- bc_ids(x) == i
            neg <- !bc_ids(x) %in% c(0, i, ms[spill_cols[[j]]])
            if (sum(pos) != 0) {
                if (method == "default") {
                    if (sum(neg) == 0) {
                        for (k in spill_cols[[j]]) {
                            spill <- 
                                mean(exprs(x)[pos, k], trim)/
                                mean(exprs(x)[pos, j], trim)
                            if (is.na(spill)) spill <- 0
                            SM[j, k] <- spill
                        }
                    } else {
                        for (k in spill_cols[[j]]) {
                            spill <- 
                                (mean(exprs(x)[pos, k], trim)- 
                                        mean(exprs(x)[neg, k], trim))/
                                (mean(exprs(x)[pos, j],  trim)-
                                        mean(exprs(x)[neg, j], trim))
                            if (is.na(spill) | spill < 0) spill <- 0
                            SM[j, k] <- spill
                        }
                    }
                } else if (method == "classic") {
                    for (k in spill_cols[[j]]) {
                        spill <- mean(exprs(x)[pos, k]/exprs(x)[pos, j], trim)
                        if (is.na(spill)) spill <- 0
                        SM[j, k] <- spill
                    }
                }
            }
        }
        colnames(SM) <- rownames(SM) <- chs
        if (interactions == "all") SM[SM < th] <- 0
        SM[bc_cols, !is.na(ms)]
    })