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
#' estimates will be set to 0.
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
#' M+/-1 (detection sensitivity), M+16 (oxide formation) and channels measuring 
#' metals that are potentially contaminated by isotopic impurites 
#' (see reference below and \code{\link{isotope_list}}).
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

#' @references 
#' Coursey, J.S., Schwab, D.J., Tsai, J.J., Dragoset, R.A. (2015).
#' Atomic weights and isotopic compositions, 
#' (available at http://physics.nist.gov/Comp).
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
        n <- ncol(es)
        chs <- colnames(es)
        
        # TODO: use helper functions to guarantee that the 
        # metals and masses are consistently parsed in the 
        # whole CATALYST package
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
            ex <- spill_cols
        } else if (interactions == "all") {
            # consider all channels
            spill_cols <- lapply(ms, function(x) which(ms != x & !is.na(ms)))
            ex <- get_spill_cols(ms, mets)
        }
        
        # compute and return compensation matrix
        SM <- matrix(diag(n), nrow=n, ncol=n, dimnames=list(chs, chs))
        for (id in ids) {
            i <- which(ms == id)
            pos <- bc_ids(x) == id
            # exclude unassigned events, i-positive & events assigned to
            # spill receiving channels from negative population
            neg <- which(!bc_ids(x) %in% c(0, id, ms[ex[[i]]]))
            pos_i <- es[pos, i]
            neg_i <- es[neg, i]
            for (j in spill_cols[[i]]) {
                pos_j <- es[pos, j]
                # further exclude events assigned to population
                # for which interaction is calculated 
                neg_j <- es[neg[bc_ids(x)[neg] != ms[j] & !(bc_ids(x)[neg] %in%  ms[ex[[j]]])], j]
                sij <- get_sij(pos_i, neg_i, pos_j, neg_j, method, trim)
                SM[i, j] <- sij
            } 
        }
        #colnames(SM) <- rownames(SM) <- chs
        SM[SM < th] <- 0
        SM[bc_cols, !is.na(ms)]
    })