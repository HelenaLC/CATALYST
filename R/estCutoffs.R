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
#' @param est     
#' logical. Should separation cutoffs be estimated? Defaults to TRUE.
#' @param verbose 
#' logical. Should extra information on progress be reported? Defaults to TRUE.
#'
#' @return
#' Will update the \code{sep_cutoffs}, \code{mhl_cutoff}, \code{counts} and 
#' \code{yields} slots of the input \code{\link{dbFrame}} and return the latter.
#' 
#' @examples
#' data(ss_beads)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_beads, y = bc_ms)
#' estCutoffs(x = re)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom stats loess

# ------------------------------------------------------------------------------

setMethod(f="estCutoffs", 
    signature=signature(x="dbFrame"), 
    definition=function(x, est=TRUE, verbose=TRUE) {
        
        seps <- seq(0, 1, .01)
        n_seps <- length(seps)
        
        ids <- rownames(x@bc_key)
        n <- length(ids)
        
        # compute well-wise yield for each separation cutoff
        if (verbose) message("Computing counts and yields...")
        yields <- counts <- matrix(0, nrow=n, ncol=n_seps)
        for (i in ids) {
            sub <- x@bc_ids == i
            ind <- ids == i
            for (j in seq_along(seps)) {
                k <- x@deltas[sub] >= seps[j]
                yields[ind, j] <- sum(k)
                counts[ind, j] <- sum(k & x@deltas[sub] < seps[j + 1])
                if (j == n_seps)
                    counts[ind, j] <- sum(k)
            }
        }
        
        # normalize each barcode to maximum yield
        norm_val <- apply(yields, 1, max)
        norm_val[norm_val == 0] <- 1
        yields <- t(sapply(1:n, function(x) yields[x, ] / norm_val[x]))
        ests <- numeric(n)
        
        # estimate cutoff parameter
        ds <- list ()
        ds[[1]] <- seps
        w <- diff(ds[[1]])[1]
        ds[[2]] <- ds[[1]][-1] - w/2
        ds[[3]] <- ds[[2]][-1] - w/2
        
        if (est) {
            if (verbose) message("Estimating separation cutoffs...")
            for (i in ids) {
                ind <- ids == i
                dy <- list()
                dy[[1]] <- yields[ind, ]
                dy[[2]] <- diff(dy[[1]]) / w
                dy[[3]] <- diff(dy[[2]]) / w
                names(dy) <- paste0("dy", 0:2)
                
                lws <- mapply( function(u, v) { 
                    predict(stats::loess(v~u, span=.3), u) }, ds, dy)
                est <- ds[[1]][which(lws[[3]] > 0 & 
                        c(lws[[3]][-1] < 0, FALSE))[1]]
                
                if (length(est) == 0 | is.na(est)) {
                    ests[ind] <- 1 
                } else {
                    ests[ind] <- est
                }
            }
        } else {
            ests <- rep(0, n)
        }
        x@sep_cutoffs <- ests
        x@counts <- counts
        x@yields <- yields
        x
    })