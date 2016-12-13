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
#'
#' @return
#' Will update the \code{bc_ids} slot of the input \code{\link{dbFrame}}
#' and return the latter.
#' 
#' @examples
#' data(ss_beads)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_beads, y = bc_ms)
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

# ------------------------------------------------------------------------------

setMethod(f="applyCutoffs", 
    signature=signature(x="dbFrame"),
    definition=function(x, mhl_cutoff=30) {
        
        # get no. of events, channel and barcode masses
        N <- nrow(x@exprs)
        nms <- colnames(x@exprs)
        ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
        bc_ms <- colnames(x@bc_key)
        
        # find which columns of loaded FCS file correspond to barcode masses
        bc_cols <- which(ms %in% bc_ms)
        n_bcs <- length(bc_cols)
        if (length(bc_cols) != n_bcs)
            stop("Not all barcode channels found.")
        
        ids <- unique(x@bc_ids)
        ids <- sort(ids[which(ids != 0)])
        
        # extract barcode columns from FCS file
        bcs <- x@exprs[, bc_cols]
        
        # compute mahalanobis distances of all events
        # given current separation cutoff
        mhl_dists <- numeric(N)
        for (i in ids) {
            inds <- which(x@bc_ids == i)
            ex <- inds[which(x@deltas[inds] < 
                    x@sep_cutoffs[which(rownames(x@bc_key) == i)])]
            inds <- inds[-which(inds %in% ex)]
            x@bc_ids[ex] <- 0
            sub  <- bcs[inds, ]
            if (length(sub) != n_bcs)
                if (nrow(sub) > n_bcs)
                    mhl_dists[inds] <- stats::mahalanobis(
                        x=sub, center=colMeans(sub), cov=stats::cov(sub))
        }
        x@bc_ids[which(mhl_dists > mhl_cutoff)] <- 0
        x@mhl_cutoff <- mhl_cutoff
        x
    })