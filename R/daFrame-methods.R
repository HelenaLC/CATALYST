# ==============================================================================
# Accessor & replacement methods for class daFrame
# ------------------------------------------------------------------------------
#' @rdname daFrame-methods
#' @title Extraction and replacement methods 
#' for objects of class \code{daFrame}
#' @aliases 
#' daFrame-methods exprs type_markers state_markers
#' n_cells cluster_codes cluster_ids
#' 
#' @description
#' Methods for accessing slots in a \code{\link{daFrame}}.
#' @return
#' \describe{
#' \item{\code{exprs}}{
#' extracts the arcsinh-transformed expressions.}
#' \item{\code{n_cells}}{
#' extracts the number of events measured per sample.}
#' \item{\code{type_markers}}{
#' extracts the antigens used for clustering.}
#' \item{\code{state_markers}}{
#' extracts antigens that were not used for clustering.}
#' \item{\code{sample_ids}}{
#' extracts the sample IDs as specified in the metadata-table.}
#' \item{\code{cluster_codes}}{
#' extracts a \code{data.frame} containing cluster codes for the 
#' \code{\link{FlowSOM}} clustering, the \code{\link{ConsensusClusterPlus}} 
#' metaclustering, and all mergings done through \code{\link{mergeClusters}}.}
#' \item{\code{cluster_ids}}{
#' extracts the numeric vector of cluster IDs 
#' as inferred by \code{\link{FlowSOM}}.}
#' }
#' 
#' @param x,object a \code{\link{daFrame}}.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
#' 
#' @examples
#' # construct daFrame
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run clustering
#' re <- cluster(re)
#' 
#' # view data summary
#' library(SummarizedExperiment)
#' cbind(metadata(re)$experiment_info, cells=n_cells(re))
#' 
#' # access row / cell data
#' head(rowData(re))
#' plot(table(cluster_ids(re)))
#' 
#' # access marker information
#' type_markers(re)
#' state_markers(re)
#' 
#' # get cluster ID correspondece between 2 clusterings
#' old_ids <- seq_len(20)
#' m <- match(old_ids, cluster_codes(re)$`meta20`)
#' new_ids <- cluster_codes(re)$`meta12`[m]
#' data.frame(old_ids, new_ids)
#' 
#' # plot relative change in area under CDF curve vs. k
#' metadata(re)$delta_area
# ------------------------------------------------------------------------------

# the following wrappers for reducedDimNames & reducedDim(s) assure that 
# reduced dimensions are called without dimnames (this is necessary b/c 
# dims. are reversed in a daFrame, and nrow(DR) != ncol(x) causes errors)
#' @rdname daFrame-methods
#' @export
setMethod("reducedDimNames", "daFrame", function(x) names(x@reducedDims))

#' @rdname daFrame-methods
#' @export
setMethod("reducedDims", "daFrame", function(x) x@reducedDims)

#' @rdname daFrame-methods
#' @param dr character string specifying the dim. reduction to extract.
#' @export
setMethod("reducedDim", "daFrame", function(x, dr = 1) 
    tryCatch(error = function(e) NULL, x@reducedDims[[dr]]))

#' @rdname daFrame-methods
#' @param value a character vector containing 
#'   the desired dimensionality reduction names.
#' @importFrom SingleCellExperiment reducedDimNames<-
#' @export
setReplaceMethod("reducedDimNames", c("daFrame", "character"),
    function(x, value) {
        dr <- reducedDims(x)
        names(dr) <- value
        names(x@int_elementMetadata) <- sprintf("idx.%s", names(dr))
        x@reducedDims <- dr
        return(x)
    }
)

#' @rdname daFrame-methods
#' @importFrom Biobase exprs
#' @export
setMethod(f="exprs",
    signature="daFrame",
    definition=function(object) return(assays(object)$exprs))

#' @rdname daFrame-methods
setMethod(f="n_cells",
    signature="daFrame",
    definition=function(x) return(metadata(x)$n_cells))

#' @rdname daFrame-methods
setMethod(f="marker_classes",
    signature="daFrame",
    definition=function(x) 
        return(setNames(unlist(colData(x)$marker_class), colnames(x))))

#' @rdname daFrame-methods
setMethod(f="type_markers",      
    signature="daFrame",
    definition=function(x) return(colnames(x)[
        colData(x)$marker_class == "type"]))

#' @rdname daFrame-methods
setMethod(f="state_markers",      
    signature="daFrame",
    definition=function(x) return(colnames(x)[
        colData(x)$marker_class == "state"]))

#' @rdname daFrame-methods
setMethod(f="sample_ids",  
    signature="daFrame", 
    definition=function(x) return(rowData(x)$sample_id))

#' @rdname daFrame-methods
setMethod(f="cluster_codes",  
    signature="daFrame", 
    definition=function(x) return(metadata(x)$cluster_codes))

#' @rdname daFrame-methods
setMethod(f="cluster_ids",  
    signature="daFrame", 
    definition=function(x) return(rowData(x)$cluster_id))