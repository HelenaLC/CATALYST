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
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
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
        x <- getMahalDists(x)
        x@bc_ids[which(x@mhl_dists > mhl_cutoff)] <- 0
        x@mhl_cutoff <- mhl_cutoff
        x
    })
