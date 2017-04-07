# ==============================================================================
# Estimate separation cutoffs
# ------------------------------------------------------------------------------

#' @rdname estCutoffs
#' @title Estimation of distance separation cutoffs
#' 
#' @description 
#' For each barcode, estimates a cutoff parameter for the 
#' distance between positive and negative barcode populations.
#'
#' @param x       
#' a \code{\link{dbFrame}}.
#' @param verbose 
#' logical. Should extra information on progress be reported? Defaults to TRUE.
#'
#' @return
#' Will update the \code{sep_cutoffs}, \code{mhl_cutoff}, \code{counts} and 
#' \code{yields} slots of the input \code{\link{dbFrame}} and return the latter.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' estCutoffs(x = re)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom stats loess

# ------------------------------------------------------------------------------

setMethod(f="estCutoffs", 
    signature=signature(x="dbFrame"), 
    definition=function(x, verbose=TRUE) {
        
        n_bcs <- nrow(bc_key(x))
        ests <- numeric(n_bcs)
        
        ds <- list ()
        ds[[1]] <- seq(0, 1, .01)
        w <- diff(ds[[1]])[1]
        ds[[2]] <- ds[[1]][-1] - w/2
        ds[[3]] <- ds[[2]][-1] - w/2
        
        if (verbose) message("Estimating separation cutoffs...")
        for (i in seq_len(n_bcs)) {
            dy <- list()
            dy[[1]] <- yields(x)[i, ]
            dy[[2]] <- diff(dy[[1]]) / w
            dy[[3]] <- diff(dy[[2]]) / w
            names(dy) <- paste0("dy", 0:2)
            
            lws <- mapply( function(u, v) { 
                predict(stats::loess(v~u, span=.3), u) }, ds, dy)
            est <- ds[[1]][which(lws[[3]] > 0 & 
                    c(lws[[3]][-1] < 0, FALSE))[1]]
            
            if (length(est) == 0 | is.na(est)) {
                ests[i] <- 1 
            } else {
                ests[i] <- est
            }
        }
        names(ests) <- rownames(bc_key(x))
        sep_cutoffs(x) <- ests
        x
    })