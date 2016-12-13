# ==============================================================================
# Compute debarcoding
# ------------------------------------------------------------------------------

#' @rdname debarcode
#' @title Compute debarcoding
#' 
#' @description 
#' For each barcode, estimates a cutoff parameter for the 
#' distance between positive and negative barcode populations.
#'
#' @param x       
#' a \code{\link{flowFrame}}.
#' @param y 
#' the debarcoding scheme. A binary matrix with sample names as row names and 
#' numeric masses as column names OR a vector of numeric masses corresponding to 
#' barcode channels. When the latter is supplied, \code{assignPrelim} will 
#' create a scheme of the appropriate format internally.
#' @param out_path
#' a character string. All outputs (population-wise FCS files and 
#' diagnostic plots) will be generated in this location. 
#' @param cofactor 
#' cofactor used for asinh transformation.
#' @param mhl_cutoff 
#' mahalanobis distance threshold above which events should be unassigned.
#' @param verbose  
#' logical. Should extra information on progress be reported? Defaults to TRUE.
#' 
#' @return
#' Will generate population-wise FCS files, and events and yields plots
#' in the specified \code{out_path}.
#' 
#' @examples
#' data(ss_beads)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_beads, y = bc_ms)
#' estCutoffs(x = re)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @export

# ------------------------------------------------------------------------------

setMethod(f="debarcode",
    signature=signature(x="flowFrame", y="data.frame"),
    definition=function(x, y, out_path, 
        cofactor=10, mhl_cutoff=30, verbose=TRUE) {
    
    # assign preliminary barcode IDs
    re <- assignPrelim(x=x, y=y, cofactor=cofactor, verbose=verbose)
    plotEvents(x=re, out_path=out_path, name_ext="_prelim")
    
    # estimate distance separation cutoffs,
    # get counts and yields
    re <- estCutoffs(x = re, verbose=verbose)
    which <- sort(unique(re@bc_ids))
    which <- which[!which == 0]
    plotYields(x=re, which_bc=c(0, which), out_path=out_path)
    
    # apply deconcolution parameters
    re <- applyCutoffs(x=re, mhl_cutoff=mhl_cutoff)
    plotEvents(x=re, out_path=out_path, name_ext="_final")

    # output population-wise FCS files
    outFCS(x=re, out_path=out_path)
    })

# --------------------------------------------------------------------------------

#' @rdname debarcode
setMethod(f="debarcode",
    signature=signature(x="flowFrame", y="vector"),
    definition=function(x, y, out_path, 
        cofactor=10, mhl_cutoff=30, verbose=TRUE) {
        n <- length(y)
        y <- data.frame(matrix(diag(n), ncol=n, 
            dimnames=list(y, y)), check.names=FALSE)
        debarcode(x, y, out_path, 
            cofactor=10, mhl_cutoff=30, verbose=TRUE)
    })

    




