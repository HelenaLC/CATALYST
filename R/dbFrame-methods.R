# ==============================================================================
# Accessor and replacement methods for class dbFrame
# ------------------------------------------------------------------------------

#' @rdname dbFrame-methods
#' @title 
#' Extraction and replacement methods for objects of class \code{dbFrame}
#' @aliases 
#' dbFrame-methods exprs bc_key bc_ids deltas normed_bcs mhl_dists 
#' sep_cutoffs sep_cutoffs<- mhl_cutoff mhl_cutoff<- counts yields
#' 
#' @description
#' Methods for replacing and accessing slots in a \code{\link{dbFrame}}.
#' @return
#' \describe{
#' \item{\code{exprs}}{extracts the raw data intensities.}
#' \item{\code{bc_key}}{extracts the barcoding scheme.}
#' \item{\code{bc_ids}}{extracts currently made event assignments.}
#' \item{\code{deltas}}{extracts barcode separations computed from normalized 
#'                      intensities. \code{sep_cutoffs} apply to these values 
#'                      (see \code{\link{applyCutoffs}}).}
#' \item{\code{normed_bcs}}{extracts normalized barcode intensities 
#'                          (see \code{\link{assignPrelim}}).}
#' \item{\code{sep_cutoffs}, \code{sep_cutoffs<-}}{extracts or replaces 
#' separation cutoffs. If option \code{sep_cutoffs} is not specified, these will
#' be used by \code{\link{applyCutoffs}}. Replacement value must be a non-
#' negative numeric with length one or same length as the number of barcodes.}
#' \item{\code{mhl_cutoff}, \code{mhl_cutoff<-}}{extracts or replaces the 
#' Mahalanobis distance threshold above which events are to be unassigned.
#' Replacement value must be a single non-negative and non-zero numeric.}
#' \item{\code{counts}}{extract the counts matrix (see \code{\link{dbFrame}}).}
#' \item{\code{yields}}{extract the yields matrix (see \code{\link{dbFrame}}).}
#' }
#' @param x,object a \code{\link{dbFrame}}.
#' @param value the replacement value.
#' 
#' @examples 
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' 
#' # set global cutoff parameter
#' sep_cutoffs(re) <- 0.4
#' re <- applyCutoffs(x = re)
#' 
#' # subset a specific population, e.g. A1: 111000
#' a1 <- bc_ids(re) == "A1"
#' head(exprs(sample_ff[a1, ]))
#' 
#' # subset unassigned events
#' unassigned <- bc_ids(re) == 0
#' head(exprs(sample_ff[unassigned, ]))
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}

# ==============================================================================
# Access slots in a dbFrame
# ------------------------------------------------------------------------------
#' 
#' @rdname dbFrame-methods
#' @importFrom flowCore exprs
#' @export
setMethod(f="exprs",
    signature="dbFrame",
    definition=function(object) return(object@exprs))

#' @rdname dbFrame-methods
setMethod(f="bc_key",      
    signature="dbFrame", 
    definition=function(x) return(x@bc_key))

#' @rdname dbFrame-methods
setMethod(f="bc_ids",      
    signature="dbFrame", 
    definition=function(x) return(x@bc_ids)) 

#' @rdname dbFrame-methods
setMethod(f="deltas",      
    signature="dbFrame", 
    definition=function(x) return(x@deltas))

#' @rdname dbFrame-methods
setMethod(f="normed_bcs",  
    signature="dbFrame", 
    definition=function(x) return(x@normed_bcs))

#' @rdname dbFrame-methods
setMethod(f="mhl_dists",  
    signature="dbFrame", 
    definition=function(x) return(x@mhl_dists))

#' @rdname dbFrame-methods
setMethod(f="sep_cutoffs", 
    signature="dbFrame", 
    definition=function(x) return(x@sep_cutoffs))

#' @rdname dbFrame-methods
setMethod(f="mhl_cutoff",  
    signature="dbFrame", 
    definition=function(x) return(x@mhl_cutoff))

#' @rdname dbFrame-methods
#' @export
setMethod(f="counts",
    signature="dbFrame",
    definition=function(x) return(x@counts))

#' @rdname dbFrame-methods
setMethod(f="yields",   
    signature="dbFrame", 
    definition=function(x) return(x@yields))

# ==============================================================================
# Replace method for slot 'bc_ids' (only used internally)
# ------------------------------------------------------------------------------

setReplaceMethod(f="bc_ids", 
    signature=signature(x="dbFrame"), 
    definition=function(x, value) {
        valid_ids <- c(0, rownames(bc_key(x)))
        if (!any(value %in% valid_ids)) {
            invalid <- value[!value %in% valid_ids]
            if (length(invalid) == 1) 
                stop("\n", invalid, " is not a valid barcode ID.",
                    "\n'bc_ids' should be either 0 = \"unassigned\"",
                    "\nor occur as rownames in the 'bc_key'.")
            if (length(invalid) > 1) 
                stop("\nBarcode IDs ", paste0(invalid, collapse=", "), 
                    " are invalid.\n'bc_ids' should be either 0 = \"",
                    "unassigned\"\nor occur as rownames in the 'bc_key'.")
        }
        x@bc_ids <- value
        return(x)
    })

# ==============================================================================
# Replace method for slot 'mhl_dists' (only used internally)
# ------------------------------------------------------------------------------

setReplaceMethod(f="mhl_dists", 
    signature=signature(x="dbFrame", value="numeric"), 
    definition=function(x, value) {
        x@mhl_dists <- value
        return(x)
    })

# ==============================================================================
# Replace method for slot 'mhl_cutoff'
# ------------------------------------------------------------------------------

#' @rdname dbFrame-methods
#' @export
setReplaceMethod(f="mhl_cutoff", 
    signature=signature(x="dbFrame", value="numeric"), 
    definition=function(x, value) {
        if (length(value) != 1)
            stop("Replacement value must be of length one.")
        if (any(value < 0)) 
            stop("Replacement value must be non-negative.")
        if (value == 0) 
            stop("Applying this cutoff will have all events unassigned.") 
        x@mhl_cutoff <- value
        return(x)
    })

#' @rdname dbFrame-methods
#' @export
setReplaceMethod(f="mhl_cutoff", 
    signature=signature(x="dbFrame", value="ANY"), 
    definition=function(x, value) {
        stop("Replacement value must be a non-negative numeric of length one.") 
    })

# ==============================================================================
# Replace method for slot 'sep_cutoffs'
# ------------------------------------------------------------------------------

#' @rdname dbFrame-methods
#' @export
setReplaceMethod(f="sep_cutoffs", 
    signature=signature(x="dbFrame", value="numeric"), 
    definition=function(x, value) {
        if (any(value < 0))
            stop("Replacement value(s) must be non-negative.")
        if (length(value) == 1) {
            x@sep_cutoffs <- rep(value, nrow(bc_key(x)))
            return(x)
        }
        if (length(value) != nrow(x@bc_key))
            stop("'Replacement value' must be of length one\n or same length",
                " as the number of rows in the 'bc_key'.")
        x@sep_cutoffs <- value
        names(x@sep_cutoffs) <- rownames(bc_key(x))
        return(x)
    })

#' @rdname dbFrame-methods
#' @export
setReplaceMethod(f="sep_cutoffs", 
    signature=signature(x="dbFrame", value="ANY"), 
    definition=function(x, value) {
        stop("Replacement value must be a non-negative numeric with length one",
            "\n or same length as the number of rows in the 'bc_key'.")
    })
