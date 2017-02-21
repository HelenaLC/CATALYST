# ==============================================================================
# Apply separation and mahalanobies distance cutoffs
# ------------------------------------------------------------------------------

#' @rdname applyCutoffs
#' @title Single-cell debarcoding (2)
#' 
#' @description Applies separation and mahalanobies distance cutoffs.
#'
#' @param x 
#' a \code{\link{dbFrame}}.
#' @param mhl_cutoff 
#' mahalanobis distance threshold above which events should be unassigned.
#' @param sep_cutoffs
#' non-negative numeric of length one or same length as the number of rows in
#' the \code{bc_key}. Specifies the distance separation cutoffs between positive 
#' and negative barcode populations above which events should be unassigned.
#' If \code{NULL} (default), \code{applyCutoffs} will try to access the 
#' 'sep_cutoffs' slot of the supplied \code{\link{dbFrame}}.
#'
#' @return 
#' Will update the \code{bc_ids} and, if specified, 
#' \code{sep_cutoffs} of the input \code{\link{dbFrame}}.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' 
#' # use global separation cutoff
#' applyCutoffs(x = re, sep_cutoffs = 0.2)
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
        if (!is.null(sep_cutoffs)) {
            sep_cutoffs(x) <- sep_cutoffs
        } else {
            if (length(sep_cutoffs(x)) == 0) 
                stop("'sep_cutoffs' need to be supplied.\n Please run",
                    " 'estCutoffs()' first, or specify cutoffs manually.")
        }
        
        # get channel and barcode masses
        nms <- colnames(x@exprs)
        ms <- as.numeric(gsub("[[:alpha:][:punct:]]", "", nms))
        
        # find which columns correspond to barcode masses
        bc_cols <- which(ms %in% colnames(x@bc_key))
        n_bcs <- length(bc_cols)
        
        # extract barcode columns from FCS file
        bcs <- x@exprs[, bc_cols]
        
        ids <- unique(x@bc_ids)
        ids <- ids[which(ids != 0)]

        # compute mahalanobis distances given current separation cutoff
        mhl_dists <- numeric(nrow(x@exprs))
        for (i in ids) {
            inds <- which(x@bc_ids == i)
            ex <- inds[x@deltas[inds] < x@sep_cutoffs[rownames(x@bc_key) == i]]
            inds <- inds[!(inds %in% ex)]
            x@bc_ids[ex] <- 0
            sub  <- bcs[inds, ]
            if (length(sub) != n_bcs)
                if (nrow(sub) > n_bcs)
                    mhl_dists[inds] <- stats::mahalanobis(
                        x=sub, center=colMeans(sub), cov=stats::cov(sub))
        }
        bc_ids(x)[mhl_dists > mhl_cutoff] <- 0
        x@mhl_dists <- mhl_dists
        mhl_cutoff(x) <- mhl_cutoff
        x
    })