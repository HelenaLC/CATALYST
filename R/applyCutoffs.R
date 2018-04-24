#' @rdname applyCutoffs
#' @title Single-cell debarcoding (2)
#' 
#' @description Applies separation and mahalanobies distance cutoffs.
#'
#' @param x 
#'   a \code{\link{dbFrame}}.
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
#' @return 
#' Will update the \code{bc_ids} and, if not already specified, 
#' \code{sep_cutoffs} & \code{mhl_cutoff} slots of \code{x}.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
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
#' @importFrom stats cov mahalanobis
# ------------------------------------------------------------------------------

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
        
        # get indices for each barcode population
        ids <- unique(bc_ids(x))
        ids <- ids[ids != 0]
        inds <- setNames(lapply(ids, function(id)
            which(bc_ids(x) == id)), ids)
        
        # get indices of events that fall below 
        # sep_cutoff & assign bc_id 0 to these
        below_th <- setNames(lapply(ids, function(id) {
            th <- sep_cutoffs(x)[rownames(bc_key(x)) == id]
            deltas(x)[inds[[id]]] < th
        }), ids)
        inds <- setNames(lapply(ids, function(id) 
            inds[[id]][!below_th[[id]]]), ids)
        for (id in ids) bc_ids(x)[inds[[id]][below_th[[id]]]] <- 0
        
        # subset barcode populations
        ms <- gsub("[[:alpha:][:punct:]]", "", colnames(exprs(x)))
        bc_cols <- which(ms %in% colnames(bc_key(x)))
        n_bcs <- length(bc_cols)
        bcs <- exprs(x)[, bc_cols]
        bcs <- setNames(lapply(inds, function(x) bcs[x, ]), names(inds))
        
        # compute mhl_dists
        mhl_dists <- numeric(nrow(exprs(x)))
        for (id in ids) {
            sub <- bcs[[id]]
            test <- (length(sub) != n_bcs) && (nrow(sub) > n_bcs)
            if (test) {
                covMat <- stats::cov(sub)
                # check if covariance matrix is invertible
                test <- tryCatch(
                    solve(covMat) %*% covMat, 
                    error=function(e) e)
                if (!inherits(test, "error"))
                    mhl_dists[inds[[id]]] <- stats::mahalanobis(
                        x=sub, center=colMeans(sub), cov=covMat)
            }
        }
        bc_ids(x)[mhl_dists > mhl_cutoff] <- 0
        mhl_dists(x) <- mhl_dists
        mhl_cutoff(x) <- mhl_cutoff
        return(x)
    }
)