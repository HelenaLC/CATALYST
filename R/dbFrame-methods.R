# ================================================================================
# class to access slots in a dbFrame
# --------------------------------------------------------------------------------

#' @rdname dbFrame-class
#' @import flowCore
#' @importFrom flowCore exprs
#' @aliases exprs exprs-method exprs,dbFrame-method dbFrame-methods
setMethod(f="exprs",       
          signature="dbFrame", 
          definition=function(object) return(object@exprs))

#' @rdname dbFrame-class
#' @aliases bc_key bc_key-method bc_key,dbFrame-method dbFrame-methods
setMethod(f="bc_key",      
          signature="dbFrame", 
          definition=function(object) return(object@bc_key))

#' @rdname dbFrame-class
#' @aliases bc_ids bc_ids-method bc_ids,dbFrame-method dbFrame-methods
setMethod(f="bc_ids",      
          signature="dbFrame", 
          definition=function(object) return(object@bc_ids))

#' @rdname dbFrame-class
#' @aliases deltas deltas-method deltas,dbFrame-method dbFrame-methods
setMethod(f="deltas",      
          signature="dbFrame", 
          definition=function(object) return(object@deltas))

#' @rdname dbFrame-class
#' @aliases normed_bcs normed_bcs-method normed_bcs,dbFrame-method dbFrame-methods
setMethod(f="normed_bcs",  
          signature="dbFrame", 
          definition=function(object) return(object@normed_bcs))

#' @rdname dbFrame-class
#' @aliases sep_cutoffs sep_cutoffs-method sep_cutoffs,dbFrame-method dbFrame-methods
setMethod(f="sep_cutoffs", 
          signature="dbFrame", 
          definition=function(object) return(object@sep_cutoffs))

#' @rdname dbFrame-class
#' @aliases mhl_cutoff mhl_cutoff-method mhl_cutoff,dbFrame-method dbFrame-methods
setMethod(f="mhl_cutoff",  
          signature="dbFrame", 
          definition=function(object) return(object@mhl_cutoff))

#' @rdname dbFrame-class
#' @aliases counts counts-method counts,dbFrame-method dbFrame-methods
setMethod(f="counts",      
          signature="dbFrame", 
          definition=function(object) return(object@counts))

#' @rdname dbFrame-class
#' @aliases yields yields-method yields,dbFrame-method dbFrame-methods
setMethod(f="yields",      
          signature="dbFrame", 
          definition=function(object) return(object@yields))

# ================================================================================
# Replace method for slot bc_ids
# --------------------------------------------------------------------------------
setReplaceMethod(f="bc_ids", 
                 signature=signature(object="dbFrame", value="numeric"), 
                 definition=function(object, value) {
                     if (!is.numeric(value))
                         stop("Replacement value for the 'bc_ids' slot of a 'dbFrame'\n",
                              "  must be a numeric.")
                     
                     if (length(value) != length(object@bc_ids))
                         stop("Replacement value for the 'bc_ids' slot of a 'dbFrame'\n",
                              "  must be of same length.")
                     
                     if (any(!(value %in% c(0, as.numeric(rownames(object@bc_ids))))))
                         stop("Replacement value for the 'bc_ids' slot of a 'dbFrame'\n",
                              "  must be in agreement with the 'bc_key'.")
                     
                     object@bc_ids <- value
                     re <- estCutoffs(x = re)
                     
                     return(object)
                 })

# throw informative error when trying to replace with anything other than numeric
setReplaceMethod(f="bc_ids",
                 signature=signature(object="dbFrame", value="ANY"),
                 definition=function(object, value)
                     stop("Replacement value for the 'bc_ids' slot of a 'dbFrame'\n",
                          "  must be a numeric of same length and with values occuring\n",
                          "  as row names in the 'bc_key'.")
)

# ================================================================================
# Replace method for slot sep_cutoffs
# --------------------------------------------------------------------------------
setReplaceMethod(f="sep_cutoffs", 
                 signature=signature(object="dbFrame", value="numeric"), 
                 definition=function(object, value) {
                     if (!is.numeric(value))
                         stop("Replacement value for the 'sep_cutoffs' slot of a 'dbFrame'\n",
                              "  must be a numeric.")
                     
                     if (length(value) != nrow(object@bc_key))
                         stop("Replacement value for the 'sep_cutoffs' slot of a 'dbFrame'\n",
                              "  must be of same length as the number of rows in slot 'bc_key'.")
                     
                     if (any(value < 0))
                         stop("Replacement value for the 'sep_cutoffs' slot of a 'dbFrame'\n",
                              "  must be non-negative.")
                     
                     object@sep_cutoffs <- value
                     object <- applyCutoffs(x = object)
                     return(object)
                 })





