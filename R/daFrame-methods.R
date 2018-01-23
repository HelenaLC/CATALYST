# ==============================================================================
# Accessor & replacement methods for class daFrame
# ------------------------------------------------------------------------------
#' @rdname daFrame-methods
#' @title 
#' Extraction and replacement methods for objects of class \code{daFrame}
#' @aliases 
#' daFrame-methods exprs lineage functional conditions 
#' sample_ids cluster_ids merging_codes merging_ids
#' 
#' @description
#' Methods for replacing and accessing slots in a \code{\link{daFrame}}.
#' @return
#' 
#' @param x,object a \code{\link{daFrame}}.
#' 
#' @examples 
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}

# ==============================================================================
# Access slots in a daFrame
# ------------------------------------------------------------------------------
#' @rdname daFrame-methods
#' @importFrom Biobase exprs
#' @export
setMethod(f="exprs",
    signature="daFrame",
    definition=function(object) return(unlist(assays(object))))

#' @rdname daFrame-methods
setMethod(f="lineage",      
    signature="daFrame", 
    definition=function(x) return(colnames(x)[as.logical(x$lineage)]))

#' @rdname daFrame-methods
setMethod(f="functional",      
    signature="daFrame",
    definition=function(x) return(colnames(x)[as.logical(x$functional)])) 

#' @rdname daFrame-methods
#' @importFrom BiocGenerics conditions
setMethod(f="conditions",      
    signature="daFrame", 
    definition=function(x) return(rowData(x)$conditions))

#' @rdname daFrame-methods
setMethod(f="sample_ids",  
    signature="daFrame", 
    definition=function(x) return(rowData(x)$sample_id))

#' @rdname daFrame-methods
setMethod(f="cluster_ids",  
    signature="daFrame", 
    definition=function(x) return(rowData(x)$cluster_id))

#' @rdname daFrame-methods
setMethod(f="cluster_codes",  
    signature="daFrame", 
    definition=function(x) return(metadata(x)$cluster_codes))