# ==============================================================================
# Apply separation and mahalanobies distance cutoffs
# ------------------------------------------------------------------------------

#' @rdname applyCutoffs
#' @title Single-cell debarcoding (2)
#' 
#' @description Applies separation and mahalanobies distance cutoffs.
#'
#' @param x a \code{\link{dbFrame}}.
#' @param mhl_cutoff mahalanobis distance threshold above which events should 
#' be unassigned. This argument will be ignored if the \code{mhl_cutoff} slot 
#' of the input \code{dbFrame} is specified.
#' @param sep_cutoffs non-negative numeric of length one or same length as the 
#' number of rows in the \code{bc_key}. Specifies the distance separation 
#' cutoffs between positive and negative barcode populations above which events 
#' should be unassigned. If \code{NULL} (default), \code{applyCutoffs} will try 
#' to access the 'sep_cutoffs' slot of the supplied \code{dbFrame}.
#'
#' @return 
#' Will update the \code{bc_ids} and, if not already specified, 
#' \code{sep_cutoffs} and \code{mhl_cutoff} slots of the input \code{dbFrame}.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' 
#' # use global separation cutoff
#' applyCutoffs(x = re, sep_cutoffs = 0.4)
#' 
#' # estimate population-specific cutoffs
#' re <- estCutoffs(x = re)
#' applyCutoffs(x = re)
#'
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom stats cov mahalanobis
# ==============================================================================

setMethod(f="applyCutoffs", 
    signature=signature(x="dbFrame"),
    definition=function(x, mhl_cutoff=30, sep_cutoffs=NULL) {
        
        # if specified, access 'mhl_cutoff' and 'sep_cutoffs' slots
        if (length(mhl_cutoff(x)) != 0) {
            mhl_cutoff <- mhl_cutoff(x)
        } else {
            # check validity of input 'mhl_cutoff'
            if (!is.numeric(mhl_cutoff) | length(mhl_cutoff) != 1)
                stop("'mhl_cutoff' must be a single ",
                    "non-negative and non-zero numeric.")
        }
        if (!is.null(sep_cutoffs)) {
            sep_cutoffs(x) <- sep_cutoffs
        } else {
            if (length(sep_cutoffs(x)) == 0) 
                stop("'sep_cutoffs' need to be supplied.\n Please run",
                    " 'estCutoffs()' first, or specify cutoffs manually.")
        }
        
        # find which columns correspond to barcode masses
        # and extract barcode columns
        ms <- gsub("[[:alpha:][:punct:]]", "", colnames(exprs(x)))
        bc_cols <- which(ms %in% colnames(bc_key(x)))
        n_bcs <- length(bc_cols)
        bcs <- exprs(x)[, bc_cols]
        
        ids <- unique(bc_ids(x))
        ids <- ids[ids != 0]
        
        # compute mahalanobis distances given current separation cutoff
        mhl_dists <- numeric(nrow(exprs(x)))
        for (i in ids) {
            inds <- which(bc_ids(x) == i)
            ex <- inds[deltas(x)[inds] < 
                    sep_cutoffs(x)[rownames(bc_key(x)) == i]]
            inds <- inds[!(inds %in% ex)]
            bc_ids(x)[ex] <- 0
            sub  <- bcs[inds, ]
            test <- (length(sub) != n_bcs) && (nrow(sub) > n_bcs)
            if (test) {
                covMat <- stats::cov(sub)
                # check if covariance matrix is invertible
                test <- tryCatch(
                    solve(covMat) %*% covMat, 
                    error=function(e) e)
                if (!inherits(test, "error"))
                    mhl_dists[inds] <- stats::mahalanobis(
                        x=sub, center=colMeans(sub), cov=covMat)
            }
        }
        bc_ids(x)[mhl_dists > mhl_cutoff] <- 0
        mhl_dists(x)  <- mhl_dists
        mhl_cutoff(x) <- mhl_cutoff
        x
    })