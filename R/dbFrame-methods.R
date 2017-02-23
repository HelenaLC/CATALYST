# ==============================================================================
# Access slots in a dbFrame
# ------------------------------------------------------------------------------

#' @rdname dbFrame-methods
#' @title dbFrame-methods
#' @importFrom flowCore exprs
#' @aliases exprs exprs-method dbFrame-method dbFrame-methods
setMethod(f="exprs",       
    signature="dbFrame", 
    definition=function(object) return(object@exprs))

#' @rdname dbFrame-methods
#' @aliases bc_key bc_key-method dbFrame-method dbFrame-methods
setMethod(f="bc_key",      
    signature="dbFrame", 
    definition=function(object) return(object@bc_key))

#' @rdname dbFrame-methods
#' @aliases bc_ids bc_ids-method dbFrame-method dbFrame-methods
setMethod(f="bc_ids",      
    signature="dbFrame", 
    definition=function(object) return(object@bc_ids))

#' @rdname dbFrame-methods
#' @aliases deltas deltas-method dbFrame-method dbFrame-methods
setMethod(f="deltas",      
    signature="dbFrame", 
    definition=function(object) return(object@deltas))

#' @rdname dbFrame-methods
#' @aliases normed_bcs normed_bcs-method dbFrame-method dbFrame-methods
setMethod(f="normed_bcs",  
    signature="dbFrame", 
    definition=function(object) return(object@normed_bcs))

#' @rdname dbFrame-methods
#' @aliases mhl_dists mhl_dists-method dbFrame-method dbFrame-methods
setMethod(f="mhl_dists",  
    signature="dbFrame", 
    definition=function(object) return(object@mhl_dists))

#' @rdname dbFrame-methods
#' @aliases sep_cutoffs sep_cutoffs-method dbFrame-method dbFrame-methods
setMethod(f="sep_cutoffs", 
    signature="dbFrame", 
    definition=function(object) return(object@sep_cutoffs))

#' @rdname dbFrame-methods
#' @aliases mhl_cutoff mhl_cutoff-method dbFrame-method dbFrame-methods
setMethod(f="mhl_cutoff",  
    signature="dbFrame", 
    definition=function(object) return(object@mhl_cutoff))

#' @rdname dbFrame-methods
#' @aliases counts counts-method dbFrame-method dbFrame-methods
setMethod(f="counts", 
    signature="dbFrame", 
    definition=function(object) return(object@counts))

#' @rdname dbFrame-methods
#' @aliases yields yields-method dbFrame-method dbFrame-methods
setMethod(f="yields",   
    signature="dbFrame", 
    definition=function(object) return(object@yields))

# ==============================================================================
# Replace method for slot 'bc_ids' (only used internally)
# ------------------------------------------------------------------------------

setReplaceMethod(f="bc_ids", 
    signature=signature(object="dbFrame"), 
    definition=function(object, value) {
        object@bc_ids <- value
        return(object)
    })

# ==============================================================================
# Replace method for slot 'mhl_dists' (only used internally)
# ------------------------------------------------------------------------------

setReplaceMethod(f="mhl_dists", 
    signature=signature(object="dbFrame", value="numeric"), 
    definition=function(object, value) {
        object@mhl_dists <- value
        return(object)
    })

# ==============================================================================
# Replace method for slot 'mhl_cutoff'
# ------------------------------------------------------------------------------

#' @rdname dbFrame-methods
#' @aliases <- dbFrame-method dbFrame-methods
#' @docType methods
#' @export
setReplaceMethod(f="mhl_cutoff", 
    signature=signature(object="dbFrame", value="numeric"), 
    definition=function(object, value) {
        if (length(value) != 1)
            stop("Replacement value must be of length one.")
        if (any(value < 0)) 
            stop("Replacement value must be non-negative.")
        if (value == 0) 
            stop("Applying this cutoff will have all events unassigned.") 
        object@mhl_cutoff <- value
        return(object)
    })

#' @rdname dbFrame-methods
#' @aliases <- dbFrame-method dbFrame-methods
#' @docType methods
#' @export
setReplaceMethod(f="mhl_cutoff", 
    signature=signature(object="dbFrame", value="ANY"), 
    definition=function(object, value) {
        stop("Replacement value must be a non-negative numeric of length one.") 
    })

# ==============================================================================
# Replace method for slot 'sep_cutoffs'
# ------------------------------------------------------------------------------

#' @rdname dbFrame-methods
#' @aliases <- dbFrame-method dbFrame-methods
#' @docType methods
#' @export
setReplaceMethod(f="sep_cutoffs", 
    signature=signature(object="dbFrame", value="numeric"), 
    definition=function(object, value) {
        if (any(value < 0))
            stop("Replacement value(s) must be non-negative.")
        if (length(value) == 1) {
            object@sep_cutoffs <- rep(value, nrow(object@bc_key))
            return(object)
        }
        if (length(value) != nrow(object@bc_key))
            stop("'Replacement value' must be of length one\n or same length",
                " as the number of rows in the 'bc_key'.")
        object@sep_cutoffs <- value
        return(object)
    })
                    
#' @rdname dbFrame-methods
#' @aliases <- dbFrame-method dbFrame-methods
#' @docType methods
#' @export
setReplaceMethod(f="sep_cutoffs", 
    signature=signature(object="dbFrame", value="ANY"), 
    definition=function(object, value) {
        stop("Replacement value must be a non-negative numeric with length one",
            "\n or same length as the number of rows in the 'bc_key'.")
    })