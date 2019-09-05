# ==============================================================================
# Accessor & replacement methods for class daFrame
# ------------------------------------------------------------------------------
#' @rdname SCE-utils
#' @title \code{\link{SingleCellExperiment}} convencience functions
#' @aliases exprs marker_classes type_markers state_markers
#'   n_cells sample_ids cluster_ids ei cluster_codes
#' 
#' @description Various wrappers to conviniently access slots
#'   in a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   created with \code{\link{prepData}}, and that are used 
#'   frequently during differential analysis.
#' 
#' @return
#' \describe{
#' \item{\code{ei}}{
#' extracts the experimental design table.}
#' \item{\code{n_cells}}{
#' extracts the number of events measured per sample.}
#' \item{\code{marker_classes}}{
#' extracts marker class assignments ("type", "state", "none").}
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
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k character string specifying the clustering to extract.
#'   Valid values are \code{names(cluster_codes(x))}.
#' 
#' @author Helena L Crowell
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' # view experimental design table
#' ei(sce)
#' 
#' # quick-access sample & cluster assignments
#' plot(table(sample_ids(sce)))
#' plot(table(cluster_ids(sce)))
#' 
#' # access marker information
#' marker_classes(sce)
#' type_markers(sce)
#' state_markers(sce)
#' 
#' # get cluster ID correspondece between 2 clusterings
#' old_ids <- seq_len(20)
#' m <- match(old_ids, cluster_codes(sce)$`meta20`)
#' new_ids <- cluster_codes(sce)$`meta12`[m]
#' data.frame(old_ids, new_ids)
#' 
#' # plot relative change in area under CDF curve vs. k
#' library(SingleCellExperiment)
#' metadata(sce)$delta_area

#' @rdname SCE-utils
#' @importFrom S4Vectors metadata
setMethod("ei", "SingleCellExperiment",
    function(x) {
        stopifnot("experiment_info" %in% names(metadata(x)))
        metadata(x)$experiment_info
    }
)

#' @rdname SCE-utils
setMethod("n_cells", "SingleCellExperiment",
    function(x) {
        stopifnot(c("n_cells", "sample_id") %in% names(ei(x)))
        setNames(ei(x)$n_cells, ei(x)$sample_id)
    }
)

#' @rdname SCE-utils
setMethod("marker_classes", "SingleCellExperiment",
    function(x) setNames(unlist(rowData(x)$marker_class), rownames(x)))

#' @rdname SCE-utils
setMethod("type_markers", "SingleCellExperiment",
    function(x) rownames(x)[marker_classes(x) == "type"])

#' @rdname SCE-utils
setMethod("state_markers", "SingleCellExperiment",
    function(x) rownames(x)[marker_classes(x) == "state"])

#' @rdname SCE-utils
setMethod("cluster_codes", "SingleCellExperiment", 
    function(x) metadata(x)$cluster_codes)

#' @rdname SCE-utils
setMethod("sample_ids", "SingleCellExperiment", 
    function(x) x$sample_id)

#' @rdname SCE-utils
setMethod("cluster_ids", 
    c("SingleCellExperiment", "missing"),
    function(x, k = NULL) {
        stopifnot(!is.null(cluster_codes(x)),
            "cluster_id" %in% names(colData(x)))
        k <- names(cluster_codes(x))[1]
        i <- as.numeric(as.character(x$cluster_id))
        droplevels(cluster_codes(x)[i, k])
    })

#' @rdname SCE-utils
setMethod("cluster_ids", 
    c("SingleCellExperiment", "character"),
    function(x, k = NULL) {
        stopifnot(!is.null(cluster_codes(x)),
            "cluster_id" %in% names(colData(x)))
        k <- .check_validity_of_k(x, k)
        i <- as.numeric(as.character(x$cluster_id))
        droplevels(cluster_codes(x)[i, k])
    })