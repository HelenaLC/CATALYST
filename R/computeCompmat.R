# ================================================================================
# Compute compensation matrix
# --------------------------------------------------------------------------------

#' @rdname computeCompmat
#' @title Compensation matrix
#' 
#' @description Computes the compensation matrix.
#'
#' @param x a \code{\link{dbFrame}} containing measured intensities, 
#' a numeric vector of column indices corresponding to barcode 
#' channels as debarcoding key, and barcode IDs.
#' @param method function to be used for computing spillover estimates. Defaults to \code{median}.
#'
#' @return 
#' Returns a square compensation matrix with dimensions and dimension names matching 
#' those of the input flowFrame. Only specified \code{bc_chs} are taken into consideration 
#' for computation. Spillover is assumed to be linear and is thence computed as the ratio 
#' of a positive barcode population's median intensity in affected and spilling channel, 
#' corrected for their median negative signals. Furthermore, spillover values are computed 
#' independently for each interacting pair of channels, that is, they are additive.
#' The current framework considers only potential (not all possible) interactions, that is,
#' +/-1M (detection sensitivity), same metals (isotopic impurites) and -16M (oxide formation).
#' By default, diagonal entries are set to 1.
#' 
#' @examples
#' data(ss_beads)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_beads, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' head(computeCompmat(x = re))
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom stats median

# --------------------------------------------------------------------------------

setMethod(f="computeCompmat", 
    signature=signature(x="dbFrame"), 
    definition=function(x, method=median) {
        
        ##### NEED BETTER WARNING MESSAGE HERE
        if (sum(rowSums(x@bc_key) == 1) != ncol(x@bc_key)) 
            stop("Cannot compute compensation matrix 
                from non single-staining experiment.")
        
        # get no. of channels, masses and metals
        chs <- colnames(x@exprs)
        n_chs <- length(chs)
        tmp <- gregexpr("[0-9]+", chs)
        ms <- as.numeric(regmatches(chs, tmp))
        mets <- gsub("[[:digit:]]+Di", "", chs)
        
        # get barcode IDs and barcodes masses
        ids <- unique(x@bc_ids)
        ids <- sort(ids[which(ids != 0)])
        bc_ms <- as.numeric(rownames(x@bc_key))
        
        # find which columns of loaded FCS file 
        # correspond to masses listed in barcode key
        bc_cols <- sapply(bc_ms, function(x) which(ms %in% x))
        
        # for each barcode channel, get spillover candidate channels
        # (+/-1M, -16M and channels measuring isotopes)
        spill_cols <- list()
        for (i in ids) {
            j <- bc_cols[ids == i]
            p1 <- m1 <- ox <- iso <- NULL
            if ((ms[j] + 1)  %in% ms) p1 <- chs[which(ms == (ms[j] + 1))]  
            if ((ms[j] - 1)  %in% ms) m1 <- chs[which(ms == (ms[j] - 1))]  
            if ((ms[j] + 16) %in% ms) ox <- chs[which(ms == (ms[j] + 16))] 
            if (mets[j] %in% mets[-j]) {
                iso <- chs[which(mets == mets[j])]
                iso <- iso[-which(iso %in% c(p1, m1) | iso == chs[j])]
            } 
            spill_cols[[j]] <- which(chs %in% c(m1, p1, iso, ox))
        }
        
        # compute and return compensation matrix
        SM <- matrix(0, nrow=n_chs, ncol=n_chs)
        for (i in ids) {
            j <- bc_cols[ids == i] 
            pos <- x@bc_ids == i
            neg <- !x@bc_ids %in% c(0, i, spill_cols[[j]])
            if (length(pos) & length(neg) != 0) {
                for (k in spill_cols[[j]]) {
                    spill <- 
                        (method(x@exprs[pos, k]) - method(x@exprs[neg, k])) / 
                        (method(x@exprs[pos, j]) - method(x@exprs[neg, j]))
                    if (is.na(spill)) spill <- 0
                    SM[j,k] <- spill
                }
            }
        }
        diag(SM) <- 1
        CM <- solve(SM)
        diag(CM) <- 1
        colnames(CM) <- rownames(CM) <- chs
        CM
    })