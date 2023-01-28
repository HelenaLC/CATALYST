#' @name SCE-accessors
#' @rdname SCE-accessors
#' @title \code{\link{SingleCellExperiment}} accessors
#' 
#' @aliases 
#' channels marker_classes type_markers state_markers 
#' ei n_cells sample_ids cluster_ids cluster_codes delta_area
#' 
#' @description 
#' Various wrappers to conviniently access slots
#' in a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' created with \code{\link{prepData}}, and that are used 
#' frequently during differential analysis.
#' 
#' @return
#' \describe{
#' \item{\code{ei}}{
#' extracts the experimental design table.}
#' \item{\code{n_cells}}{
#' extracts the number of events measured per sample.}
#' \item{\code{channels}}{
#' extracts the original FCS file's channel names.}
#' \item{\code{marker_classes}}{
#' extracts marker class assignments ("type", "state", "none").}
#' \item{\code{type_markers}}{
#' extracts the antigens used for clustering.}
#' \item{\code{state_markers}}{
#' extracts antigens that were not used for clustering.}
#' \item{\code{sample_ids}}{
#' extracts the sample IDs as specified in the metadata-table.}
#' \item{\code{cluster_ids}}{
#' extracts the numeric vector of cluster IDs 
#' as inferred by \code{\link{FlowSOM}}.}
#' \item{\code{cluster_codes}}{
#' extracts a \code{data.frame} containing cluster codes for the 
#' \code{\link{FlowSOM}} clustering, the \code{\link{ConsensusClusterPlus}} 
#' metaclustering, and all mergings done through \code{\link{mergeClusters}}.}
#' \item{\code{delta_area}}{
#' extracts the delta area plot stored in the 
#' SCE's \code{metadata} by \code{\link{cluster}}} 
#' }
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k character string specifying the clustering to extract.
#'   Valid values are \code{names(cluster_codes(x))}.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
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
#' # access specific clustering resolution
#' table(cluster_ids(sce, k = "meta8"))
#' 
#' # access marker information
#' channels(sce)
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
#' # view delta area plot (relative change in area 
#' # under CDF curve vs. the number of clusters 'k')
#' delta_area(sce)
NULL

#' @rdname SCE-accessors
#' @importFrom S4Vectors metadata
setMethod("ei", "SingleCellExperiment",
    function(x) .get_ei(x))

#' @rdname SCE-accessors
#' @importFrom stats setNames
setMethod("n_cells", "SingleCellExperiment",
    function(x) {
        sids <- x$sample_id
        if (is.null(sids)) return(NULL)
        table(droplevels(factor(sids)))
    })

#' @rdname SCE-accessors
#' @importFrom SummarizedExperiment rowData
setMethod("channels", "SingleCellExperiment",
    function(x) {
        chs <- rowData(x)$channel_name
        if (is.null(chs)) return(NULL)
        names(chs) <- rownames(x); chs
    })

#' @rdname SCE-accessors
#' @importFrom SummarizedExperiment rowData
setMethod("marker_classes", "SingleCellExperiment",
    function(x) {
        mcs <- rowData(x)$marker_class
        if (is.null(mcs)) return(NULL)
        names(mcs) <- rownames(x); mcs
    })

#' @rdname SCE-accessors
#' @importFrom SummarizedExperiment rowData
setMethod("type_markers", "SingleCellExperiment",
    function(x) {
        mcs <- rowData(x)$marker_class
        if (is.null(mcs)) return(NULL)
        rownames(x)[mcs == "type"]
    })

#' @rdname SCE-accessors
#' @importFrom SummarizedExperiment rowData
setMethod("state_markers", "SingleCellExperiment",
    function(x) {
        mcs <- rowData(x)$marker_class
        if (is.null(mcs)) return(NULL)
        rownames(x)[mcs == "state"]
    })

#' @rdname SCE-accessors
setMethod("sample_ids", "SingleCellExperiment", 
    function(x) x$sample_id)

#' @rdname SCE-accessors
setMethod("cluster_ids", 
    c("SingleCellExperiment", "missing"),
    function(x, k = NULL) {
        stopifnot(!is.null(x$cluster_id))
        return(x$cluster_id)
    })

#' @rdname SCE-accessors
setMethod("cluster_ids", 
    c("SingleCellExperiment", "character"),
    function(x, k = NULL) {
        stopifnot(!is.null(cluster_codes(x)), !is.null(x$cluster_id))
        k <- .check_k(x, k)
        codes <- cluster_codes(x)
        m <- match(x$cluster_id, codes[, 1])
        droplevels(codes[m, k])
    })

#' @rdname SCE-accessors
#' @importFrom S4Vectors metadata
setMethod("cluster_codes", "SingleCellExperiment", 
    function(x) metadata(x)$cluster_codes)

#' @rdname SCE-accessors
#' @importFrom S4Vectors metadata
setMethod("delta_area", "SingleCellExperiment",
    function(x) metadata(x)$delta_area)