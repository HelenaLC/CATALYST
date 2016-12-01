# ================================================================================
# Assign preliminary barcode IDs
# --------------------------------------------------------------------------------

#' @rdname assignPrelim
#' @title Single-cell debarcoding (1)
#' 
#' @description 
#' Assigns a preliminary barcode ID to each event of normalized barcode intensities 
#' and computes the barcode separation.
#'
#' @param x        a \code{\link{flowFrame}}.
#' @param y        the debarcoding scheme.
#'                 A binary matrix with numeric masses as column names OR
#'                 a vector of numeric masses corresponding to barcode channels.
#'                 When the latter is supplied, \code{assignPrelim} will create
#'                 a matrix of the appropriate format internally.
#' @param cofactor cofactor used for asinh transformation.
#' @param verbose  logical. Should extra information on progress be reported? Defaults to TRUE.
#'
#' @return 
#' Returns a \code{dbFrame} containing measured intensities, the debarcoding key, 
#' a numeric verctor of barcode IDs and separations between positive and negative 
#' barcode populations, and a matrix of normalized barcode intensities. 
#' See \code{\link{dbFrame}} for more details.

#'
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @examples
#' data(ss_beads)
#' bc_ms <- c(139, 141:156, 158:176)
#' assignPrelim(x = ss_beads, y = bc_ms)
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom stats quantile
#' @importFrom flowCore colnames exprs
#' @importClassesFrom flowCore flowFrame

# --------------------------------------------------------------------------------

setMethod(f="assignPrelim",
    signature=signature(x="flowFrame", y="data.frame"),
    definition=function(x, y, cofactor = 10, verbose = TRUE) {
        
        # get masses, intensities and no. of events
        nms <- flowCore::colnames(x)
        ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
        es <- flowCore::exprs(x)
        N <- nrow(es)
        
        # get barcodes masses and check for validity of barcode channels
        ids   <- as.numeric(rownames(y))
        bc_ms <- as.numeric(colnames(y))  
        if (any(!(bc_ms %in% ms)))
            stop("Invalid barcode channel(s) specified.")
        
        # find which columns of loaded FCS file correspond to barcode masses
        bc_cols <- sapply(bc_ms, function(x) which(ms %in% x))
        if (length(bc_cols) != length(bc_ms))
            stop("Not all barcode channels found.")
        
        # extract and transform barcode columns from FCS file
        bcs <- asinh(es[, bc_cols] / cofactor)
        
        # COMPUTE DEBARCODING
        # for each ea. of barcode intensities,
        # assign barcode ID and calculate separation
        if (verbose) message("Debarcoding data...")
        bc_ids <- debarcode(bcs, bcs, y, ids, verbose)$bc_ids
        
        # NORMALIZE BY POPULATION
        # rescale transformed barcodes for ea. population
        # using preliminary assignments
        if (verbose) message("Normalizing...")
        normed_bcs <- matrix(0, nrow=N, ncol=ncol(bcs), 
            dimnames=list(NULL, colnames(bcs)))
        for (i in ids) {
            inds <- bc_ids == i
            if (length(inds) > 1) {
                pos_bcs <- bcs[inds, ids == i]
                norm_val <- stats::quantile(pos_bcs, .95)
                normed_bcs[inds, ] <- bcs[inds, ] / norm_val
            }
        }
        
        # COMPUTE DEBARCODING
        # for normalized barcode intensities
        if (verbose) message("Debarcoding normalized data...")
        re <- debarcode(normed_bcs, bcs, y, ids, verbose)
        
        return(new(Class="dbFrame", 
            exprs=es, bc_key=y, bc_ids=re$bc_ids, 
            deltas=re$deltas, normed_bcs=normed_bcs))
    })

# --------------------------------------------------------------------------------

#' @rdname assignPrelim
setMethod(f="assignPrelim",
    signature=signature(x="flowFrame", y="vector"),
    definition=function(x, y, cofactor = 10, verbose = TRUE) {
        n <- length(y)
        y <- data.frame(matrix(diag(n), ncol=n, dimnames=list(y, y)), check.names=FALSE)
        assignPrelim(x, y, cofactor, verbose)
    })