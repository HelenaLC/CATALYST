# ==============================================================================
# Assign preliminary barcode IDs
# ------------------------------------------------------------------------------

#' @rdname assignPrelim
#' @title Single-cell debarcoding (1)
#' 
#' @description 
#' Assigns a preliminary barcode ID to each event.
#'
#' @param x 
#' a \code{\link{flowFrame}} or character of an FCS file name.
#' @param y 
#' the debarcoding scheme. A binary matrix with sample names as row names and 
#' numeric masses as column names OR a vector of numeric masses corresponding
#' to barcode channels. When the latter is supplied, \code{assignPrelim} will 
#' create a scheme of the appropriate format internally.
#' @param cofactor 
#' cofactor used for asinh transformation.
#' @param verbose
#' logical. Should extra information on progress be reported? Defaults to TRUE.
#'
#' @return 
#' Returns a \code{\link{dbFrame}} containing measurement intensities,
#' the debarcoding key, a numeric verctor of barcode IDs and separations
#' between positive and negative barcode populations, and barcode intensities
#' normalized by population. 
#'
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @examples
#' data(sample_ff, sample_key)
#' assignPrelim(x = sample_ff, y = sample_key)
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import matrixStats
#' @importFrom stats quantile
#' @importFrom flowCore colnames exprs flowFrame read.FCS
# ==============================================================================

setMethod(f="assignPrelim",
    signature=signature(x="flowFrame", y="data.frame"),
    definition=function(x, y, cofactor=10, verbose=TRUE) {
        
        # get masses, intensities and no. of events
        ms <- gsub("[[:alpha:][:punct:]]", "", flowCore::colnames(x))
        es <- flowCore::exprs(x)
        
        # get barcodes masses and check for validity of barcode channels
        n_bcs <- nrow(y)
        ids   <- rownames(y)
        bc_ms <- as.numeric(colnames(y))  
        if (any(!(bc_ms %in% ms)))
            stop("Invalid barcode channel(s) specified.")
        
        # find which columns of loaded FCS file correspond to barcode masses
        bc_cols <- vapply(bc_ms, function(x) which(ms == x), numeric(1))
        if (length(bc_cols) != length(bc_ms))
            stop("Not all barcode channels found.")
        
        # extract and transform barcode columns from FCS file
        bcs <- asinh(es[, bc_cols] / cofactor)
        
        # COMPUTE DEBARCODING
        # assign barcode ID to ea. event 
        if (verbose) message("Debarcoding data...")
        bc_ids <- get_ids(bcs, y, ids, verbose)
        
        # NORMALIZE BY POPULATION
        # rescale transformed barcodes for ea. population
        # using preliminary assignments
        if (verbose) message("Normalizing...")
        normed_bcs <- matrix(0, nrow=nrow(x), ncol=ncol(bcs), 
            dimnames=list(NULL, colnames(bcs)))
        for (i in ids) {
            key_ind <- ids == i
            inds <- bc_ids == i
            if (sum(inds) > 1) {
                pos_bcs <- bcs[inds, y[key_ind, ] == 1]
                norm_val <- stats::quantile(pos_bcs, .95)
                normed_bcs[inds, ] <- bcs[inds, ] / norm_val
            }
        }
        
        # get deltas from normalized intensities 
        if (verbose) message("Computing deltas...")
        deltas <- get_deltas(normed_bcs, y, verbose)
        
        seps <- seq(0, 1, .01)
        n_seps <- length(seps)
        
        # compute well-wise yield for each separation cutoff
        if (verbose) message("Computing counts and yields...")
        yields <- counts <- matrix(0, nrow=n_bcs, ncol=n_seps)
        for (i in seq_along(ids)) {
            sub <- bc_ids == ids[i]
            for (j in seq_along(seps)) {
                k <- deltas[sub] >= seps[j]
                yields[i, j] <- sum(k)
                counts[i, j] <- sum(k & deltas[sub] < seps[j + 1])
                if (j == n_seps)
                    counts[i, j] <- sum(k)
            }
        }
        
        # normalize yields
        norm_val <- matrixStats::rowMaxs(yields)
        norm_val[norm_val == 0] <- 1
        yields <- yields / norm_val
        
        rownames(counts) <- rownames(yields) <- ids
        colnames(counts) <- colnames(yields) <- seps
        
        new(Class="dbFrame", 
            exprs=es, bc_key=y, bc_ids=bc_ids, 
            deltas=deltas, normed_bcs=normed_bcs,
            counts=counts, yields=yields)
    })

# ------------------------------------------------------------------------------

#' @rdname assignPrelim
setMethod(f="assignPrelim",
    signature=signature(x="flowFrame", y="vector"),
    definition=function(x, y, cofactor=10, verbose=TRUE) {
        n <- length(y)
        y <- data.frame(matrix(diag(n), ncol=n, 
            dimnames=list(y, y)), check.names=FALSE)
        assignPrelim(x, y, cofactor, verbose)
    })

# ------------------------------------------------------------------------------

#' @rdname assignPrelim
setMethod(f="assignPrelim",
    signature=signature(x="character", y="data.frame"),
    definition=function(x, y, cofactor=10, verbose=TRUE) {
        if (length(x) != 1) 
            stop("'x' should be a single character or flowFrame.")
        if (sum(flowCore::isFCSfile(x)) != 1) 
            stop(x, " is not a valid FCS file.")
        x <- flowCore::read.FCS(x)
        assignPrelim(x, y, cofactor, verbose)
    })

# ------------------------------------------------------------------------------

#' @rdname assignPrelim
setMethod(f="assignPrelim",
    signature=signature(x="character", y="vector"),
    definition=function(x, y, cofactor=10, verbose=TRUE) {
        if (length(x) != 1) 
            stop("'x' should be a single character or flowFrame.")
        if (sum(flowCore::isFCSfile(x)) != 1) 
            stop(x, " is not a valid FCS file.")
        n <- length(y)
        y <- data.frame(matrix(diag(n), ncol=n, 
            dimnames=list(y, y)), check.names=FALSE)
        assignPrelim(x, y, cofactor, verbose)
    })
