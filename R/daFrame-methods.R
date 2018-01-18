# ==============================================================================
# Accessor & replacement methods for class daFrame
# ------------------------------------------------------------------------------
#' @rdname daFrame-methods
#' @title 
#' Extraction and replacement methods for objects of class \code{daFrame}
#' @aliases 
#' daFrame-methods data panel metadata cluster_ids 
#' merging_ids cluster_ids<- merging_ids<-
#' 
#' @description
#' Methods for replacing and accessing slots in a \code{\link{daFrame}}.
#' @return
#' \describe{
#' \item{\code{data}}{extracts the raw data.}
#' \item{\code{panel}}{extracts the panel.}
#' \item{\code{metadata}}{extracts the metadata.}
#' \item{\code{cluster_ids}}{extracts current cluster assignments.}
#' \item{\code{merging_ids}}{extracts second layer of cluster assignments.}
#' }
#' @param x,object a \code{\link{daFrame}}.
#' @param value the replacement value.
#' 
#' @examples 
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}

# ==============================================================================
# Access slots in a daFrame
# ------------------------------------------------------------------------------
#' @rdname daFrame-methods
setMethod(f="data",
    signature="daFrame",
    definition=function(x) return(x@data))

#' @rdname daFrame-methods
setMethod(f="panel",      
    signature="daFrame", 
    definition=function(x) return(x@panel))

#' @rdname daFrame-methods
setMethod(f="metadata",      
    signature="daFrame", 
    definition=function(x) return(x@metadata)) 

#' @rdname daFrame-methods
setMethod(f="cluster_ids",      
    signature="daFrame", 
    definition=function(x) return(x@cluster_ids))

#' @rdname daFrame-methods
setMethod(f="merging_ids",  
    signature="daFrame", 
    definition=function(x) return(x@merging_ids))