# ==============================================================================
# Generics for normalization
# ------------------------------------------------------------------------------

#' @rdname concatFCS
#' @param ... optional arguments.
#' @export
setGeneric("concatFCS", 
    function(x, ...) standardGeneric("concatFCS"))

#' @rdname normCytof
#' @param ... optional arguments.
#' @export
setGeneric("normCytof", 
    function(x, y, ...) standardGeneric("normCytof"))

# ==============================================================================
# Generics to access slots in a dbFrame
# ------------------------------------------------------------------------------

#' @rdname dbFrame-methods
#' @name dbFrame-methods
#' @export
setGeneric("bc_key", 
    function(x) standardGeneric("bc_key"))

#' @rdname dbFrame-methods
#' @export
setGeneric("bc_ids",      
    function(x) standardGeneric("bc_ids"))

#' @rdname dbFrame-methods
#' @export
setGeneric("deltas",  
    function(x) standardGeneric("deltas"))

#' @rdname dbFrame-methods
#' @export
setGeneric("normed_bcs", 
    function(x) standardGeneric("normed_bcs"))

#' @rdname dbFrame-methods
#' @export
setGeneric("mhl_dists", 
    function(x) standardGeneric("mhl_dists"))

#' @rdname dbFrame-methods
#' @export
setGeneric("sep_cutoffs", 
    function(x) standardGeneric("sep_cutoffs"))

#' @rdname dbFrame-methods
#' @export
setGeneric("mhl_cutoff",  
    function(x) standardGeneric("mhl_cutoff"))

#' @rdname dbFrame-methods
#' @export
setGeneric("counts",      
    function(x) standardGeneric("counts"))

#' @rdname dbFrame-methods
#' @export
setGeneric("yields",      
    function(x) standardGeneric("yields"))

# ==============================================================================
# Generics to replace slots in a dbFrame
# ------------------------------------------------------------------------------

setGeneric("bc_ids<-",      
    function(x, value) standardGeneric("bc_ids<-"))

setGeneric("mhl_dists<-",   
    function(x, value) standardGeneric("mhl_dists<-"))

setGeneric("sep_cutoffs<-", 
    function(x, value) standardGeneric("sep_cutoffs<-"))

setGeneric("mhl_cutoff<-",  
    function(x, value) standardGeneric("mhl_cutoff<-"))

setGeneric("counts<-",      
    function(x, value) standardGeneric("counts<-"))

setGeneric("yields<-",      
    function(x, value) standardGeneric("yields<-"))

# ==============================================================================
# Generics for debarcoding
# ------------------------------------------------------------------------------

#' @rdname assignPrelim
#' @param ... optional arguments.
#' @export
setGeneric("assignPrelim", 
    function(x, y, ...) standardGeneric("assignPrelim"))

#' @rdname estCutoffs
#' @param ... optional arguments.
#' @export
setGeneric("estCutoffs",   
    function(x, ...) standardGeneric("estCutoffs"))

#' @rdname applyCutoffs
#' @param ... optional arguments.
#' @export
setGeneric("applyCutoffs", 
    function(x, ...) standardGeneric("applyCutoffs"))

#' @rdname outFrames
#' @param ... optional arguments.
#' @export
setGeneric("outFrames",      
    function(x, ...) standardGeneric("outFrames"))

#' @rdname outFCS
#' @param ... optional arguments.
#' @export
setGeneric("outFCS",      
    function(x, y, out_path=tempdir(), ...) standardGeneric("outFCS"))

# ==============================================================================
# Generics for plotting
# ------------------------------------------------------------------------------

#' @rdname plotYields
#' @param ... optional arguments.
#' @export
setGeneric("plotYields", 
    function(x, ...) standardGeneric("plotYields"))

#' @rdname plotEvents
#' @param ... optional arguments.
#' @export
setGeneric("plotEvents", 
    function(x, ...) standardGeneric("plotEvents"))

#' @rdname plotMahal
#' @param ... optional arguments.
#' @export
setGeneric("plotMahal", 
    function(x, ...) standardGeneric("plotMahal"))

# ==============================================================================
# Generics for compensation
# ------------------------------------------------------------------------------
#' @rdname computeSpillmat
#' @param ... optional arguments.
#' @export
setGeneric("computeSpillmat", 
    function(x, ...) standardGeneric("computeSpillmat"))

#' @rdname plotSpillmat
#' @param ... optional arguments.
#' @export
setGeneric("plotSpillmat", 
    function(bc_ms, SM, ...) standardGeneric("plotSpillmat"))

#' @rdname adaptSpillmat
#' @param ... optional arguments.
#' @export
setGeneric("adaptSpillmat", 
    function(input_sm, out_chs, ...) standardGeneric("adaptSpillmat"))

#' @rdname compCytof
#' @param ... optional arguments.
#' @export
setGeneric("compCytof", 
    function(x, y, ...) standardGeneric("compCytof"))

# ==============================================================================
# Generics to access slots in a daFrame
# ------------------------------------------------------------------------------

#' @rdname daFrame-class
#' @param ... optional arguments.
#' @export
setGeneric("daFrame", 
    function(x, panel, md, ...) standardGeneric("daFrame"))

#' @rdname daFrame-methods
#' @export
setGeneric("n_cells", 
    function(x) standardGeneric("n_cells"))

#' @rdname daFrame-methods
#' @export
setGeneric("marker_classes", 
    function(x) standardGeneric("marker_classes"))

#' @rdname daFrame-methods
#' @export
setGeneric("type_markers", 
    function(x) standardGeneric("type_markers"))

#' @rdname daFrame-methods
#' @export
setGeneric("state_markers",      
    function(x) standardGeneric("state_markers"))

#' @rdname daFrame-methods
#' @export
setGeneric("sample_ids", 
    function(x) standardGeneric("sample_ids"))

#' @rdname daFrame-methods
#' @export
setGeneric("cluster_codes",  
    function(x) standardGeneric("cluster_codes"))

#' @rdname daFrame-methods
#' @export
setGeneric("cluster_ids",  
    function(x) standardGeneric("cluster_ids"))

# ==============================================================================
# Generics for differential analysis
# ------------------------------------------------------------------------------

#' @rdname guessPanel
#' @param ... optional arguments.
#' @export
setGeneric("guessPanel", 
    function(x, ...) standardGeneric("guessPanel"))

#' @rdname plotCounts
#' @param ... optional arguments.
#' @export
setGeneric("plotCounts", 
    function(x, ...) standardGeneric("plotCounts"))

#' @rdname plotExprs
#' @param ... optional arguments.
#' @export
setGeneric("plotExprs", 
    function(x, ...) standardGeneric("plotExprs"))

#' @rdname plotClusterExprs
#' @param ... optional arguments.
#' @export
setGeneric("plotClusterExprs", 
    function(x, ...) standardGeneric("plotClusterExprs"))

#' @rdname plotMDS
#' @param ... optional arguments.
#' @export
setGeneric("plotMDS", 
    function(x, ...) standardGeneric("plotMDS"))

#' @rdname plotExprHeatmap
#' @param ... optional arguments.
#' @export
setGeneric("plotExprHeatmap", 
    function(x, ...) standardGeneric("plotExprHeatmap"))

#' @rdname cluster
#' @param ... optional arguments.
#' @export
setGeneric("cluster", 
    function(x, ...) standardGeneric("cluster"))

#' @rdname mergeClusters
#' @param ... optional arguments.
#' @export
setGeneric("mergeClusters", 
    function(x, k, table, id) standardGeneric("mergeClusters"))

#' @rdname plotAbundances
#' @param ... optional arguments.
#' @export
setGeneric("plotAbundances", 
    function(x, ...) standardGeneric("plotAbundances"))

#' @rdname plotClusterHeatmap
#' @param ... optional arguments.
#' @export
setGeneric("plotClusterHeatmap", 
    function(x, ...) standardGeneric("plotClusterHeatmap"))

#' @rdname plotMedExprs
#' @param ... optional arguments.
#' @export
setGeneric("plotMedExprs", 
    function(x, ...) standardGeneric("plotMedExprs"))

#' @rdname plotCodes
#' @param ... optional arguments.
#' @export
setGeneric("plotCodes", 
    function(x, ...) standardGeneric("plotCodes"))

#' @rdname plotNRS
#' @param ... optional arguments.
#' @export
setGeneric("plotNRS", 
    function(x, ...) standardGeneric("plotNRS"))

#' @rdname tSNE
#' @param ... optional arguments.
#' @export
setGeneric("tSNE", 
    function(x, ...) standardGeneric("tSNE"))

#' @rdname plotSNE
#' @param ... optional arguments.
#' @export
setGeneric("plotSNE", 
    function(x, ...) standardGeneric("plotSNE"))

#' @rdname plotDiffHeatmap
#' @param ... optional arguments.
#' @export
setGeneric("plotDiffHeatmap", 
    function(x, y, ...) standardGeneric("plotDiffHeatmap"))

#' @rdname extractClusters
#' @param ... optional arguments.
#' @export
setGeneric("extractClusters", 
    function(x, k, ...) standardGeneric("extractClusters"))