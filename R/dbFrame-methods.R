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
# Replace method for slot bc_ids
# ------------------------------------------------------------------------------
setReplaceMethod(f="bc_ids", 
                 signature=signature(object="dbFrame", value="numeric"), 
                 definition=function(object, value) {
                     if (!is.numeric(value))
                         stop("Replacement value for the 'bc_ids' slot ",
                             "of a 'dbFrame'\n  must be numeric.")
                     
                     if (length(value) != length(object@bc_ids))
                         stop("Replacement value for the 'bc_ids' slot ",
                             "of a 'dbFrame'\n  must be of same length.")
                     
                     if (any(!(value %in% 
                             c(0, as.numeric(rownames(object@bc_key))))))
                         stop("Replacement value for the 'bc_ids' slot ",
                             "of a 'dbFrame'\n  must be in agreement with ",
                             "the specified 'bc_key'. Valid IDs are:\n  ", 
                             paste(c(0, as.numeric(rownames(object@bc_key))), 
                                 fill=""))
                     
                     object@bc_ids <- value
                     object <- estCutoffs(x = object)
                     object <- applyCutoffs(x = object)
                     return(object)
                 })

# ==============================================================================
# Replace method for slot sep_cutoffs
# ------------------------------------------------------------------------------
#' @rdname dbFrame-methods
#' @aliases <- dbFrame-method dbFrame-methods
#' @docType methods
#' @export
setReplaceMethod(f="sep_cutoffs", 
                 signature=signature(object="dbFrame", value="numeric"), 
                 definition=function(object, value) {
                     if (!is.numeric(value))
                         stop("Replacement value for the 'sep_cutoffs' slot ",
                             "of a 'dbFrame'\n  must be a numeric.")
                     
                     if (length(value) != nrow(object@bc_key))
                         stop("Replacement value for the 'sep_cutoffs' slot ",
                             "of a 'dbFrame'\n  must be of same length as the ",
                             "number of rows in slot 'bc_key'.")
                     
                     if (any(value < 0))
                         stop("Replacement value for the 'sep_cutoffs' slot ",
                             "of a 'dbFrame'\n  must be non-negative.")
                     
                     object@sep_cutoffs <- value
                     return(object)
                 })