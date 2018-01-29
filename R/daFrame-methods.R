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
#' Methods for accessing slots in a \code{\link{daFrame}}.
#' @return
#' \describe{
#' \item{\code{exprs}}{extracts the arcsinh-transformed expressions.}
#' \item{\code{n_events}}{extracts the number of events measured per sample.}
#' \item{\code{type1}}{extracts the antigens used for clustering.}
#' \item{\code{type2}}{extracts antigens that were not used for clustering.}
#' \item{\code{sample_ids}}{extracts the sample IDs 
#' as specified in the metadata-table.}
#' \item{\code{cluster_ids}}{extracts the numeric vector of cluster IDs
#' as inferred by \code{\link{FlowSOM}}.}
#' \item{\code{cluster_codes}}{extracts a \code{data.frame} containing 
#' cluster codes for the \code{\link{FlowSOM}} clustering, 
#' the \code{\link{ConsensusClusterPlus}} metaclustering, 
#' and all mergings done through \code{\link{mergeClusters}}.}
#' }
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
setMethod(f="n_events",
    signature="daFrame",
    definition=function(x) return(metadata(x)$n_events))

#' @rdname daFrame-methods
setMethod(f="type1",      
    signature="daFrame", 
    definition=function(x) return(colnames(x)[as.logical(colData(x)$type1)]))

#' @rdname daFrame-methods
setMethod(f="type2",      
    signature="daFrame",
    definition=function(x) return(colnames(x)[as.logical(colData(x)$type2)]))

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