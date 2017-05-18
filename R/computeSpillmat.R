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
#' Defaults to \code{mean}.
#' @param trim
#' trim value used for estimation of spill values. Note that
#' \code{trim = 0.5} is equivalent to \code{method = "median"}.
#'
#' @return 
#' Returns a square compensation matrix with dimensions and dimension names 
#' matching those of the input flowFrame. Spillover is assumed to be linear and 
#' is thence computed as the ratio of a positive barcode population's median 
#' or (trimmed) mean intensity in affected and spilling channel. Furthermore, 
#' on the basis of their additive nature, spillover values are computed 
#' independently for each interacting pair of channels. The current framework 
#' considers only potential (not all possible) interactions, that is,
#' M+/-1 (detection sensitivity), same metals (isotopic impurites) and M+16M 
#' (oxide formation). By default, diagonal entries are set to 1.
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
    
    definition=function(x, method="mean", trim=.08) {
        
        if (method == "mean") {
            method=function(...) mean(..., trim=trim) 
        } else if (method == "median") {
            method=function(...) median(...)
        }
        
        if (sum(rowSums(bc_key(x)) == 1) != ncol(bc_key(x))) 
            stop("Cannot compute compensation matrix 
                from non single-staining experiment.")
        
        # get no. of channels, masses and metals
        chs <- colnames(exprs(x))
        n_chs <- length(chs)
        tmp <- gregexpr("[0-9]+", chs)
        ms <- as.numeric(regmatches(chs, tmp))
        mets <- gsub("[[:digit:]]+Di", "", chs)
        
        # get barcode IDs and barcodes masses
        ids <- unique(bc_ids(x))
        ids <- sort(ids[which(ids != 0)])
        bc_ms <- as.numeric(rownames(bc_key(x)))
        
        # find which columns of loaded FCS file 
        # correspond to masses listed in barcode key
        bc_cols <- vapply(bc_ms, function(x) which(ms == x), numeric(1))
        
        # for each barcode channel, get spillover candidate channels
        # (+/-1M, -16M and channels measuring isotopes)
        spill_cols <- get_spill_cols(ms, mets)
        
        # compute and return compensation matrix
        SM <- diag(n_chs)
        for (i in ids) {
            j <- bc_cols[ids == i]
            pos <- bc_ids(x) == i
            neg <- !bc_ids(x) %in% c(0, i, ms[spill_cols[[j]]])
            if (sum(pos) != 0) {
                if (sum(neg) == 0) {
                    for (k in spill_cols[[j]]) {
                        spill <-
                            method(exprs(x)[pos, k])/method(exprs(x)[pos, j])
                        if (is.na(spill)) spill <- 0
                        SM[j, k] <- spill
                    }
                } else {
                    for (k in spill_cols[[j]]) {
                        spill <-
                            (method(exprs(x)[pos, k])-method(exprs(x)[neg, k]))/
                            (method(exprs(x)[pos, j])-method(exprs(x)[neg, j]))
                        if (is.na(spill) | spill < 0) spill <- 0
                        SM[j, k] <- spill
                    }
                }
            }
        }
        colnames(SM) <- rownames(SM) <- chs
        SM[bc_cols, !is.na(ms)]
    })