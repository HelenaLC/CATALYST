# ==============================================================================
# Generics for normalization
# ------------------------------------------------------------------------------

#' @rdname normCytof
#' @param ... optional arguments.
#' @export
setGeneric(name="normCytof", 
    package="CATALYST", 
    def=function(x, y, ...) standardGeneric("normCytof"))

#' @rdname concatFCS
#' @param ... optional arguments.
#' @export
setGeneric(name="concatFCS", 
    package="CATALYST", 
    def=function(x, ...) standardGeneric("concatFCS"))

# ==============================================================================
# Generics to access slots in a dbFrame
# ------------------------------------------------------------------------------

#' @rdname dbFrame-methods
#' @name dbFrame-methods
#' @export
setGeneric(name="bc_key", 
    package="CATALYST", 
    def=function(x) standardGeneric("bc_key"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="bc_ids",      
    package="CATALYST", 
    def=function(x) standardGeneric("bc_ids"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="deltas",  
    package="CATALYST", 
    def=function(x) standardGeneric("deltas"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="normed_bcs",  
    package="CATALYST", 
    def=function(x) standardGeneric("normed_bcs"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="mhl_dists",  
    package="CATALYST", 
    def=function(x) standardGeneric("mhl_dists"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="sep_cutoffs", 
    package="CATALYST", 
    def=function(x) standardGeneric("sep_cutoffs"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="mhl_cutoff",  
    package="CATALYST", 
    def=function(x) standardGeneric("mhl_cutoff"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="counts",      
    package="CATALYST", 
    def=function(x) standardGeneric("counts"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="yields",      
    package="CATALYST", 
    def=function(x) standardGeneric("yields"))

# ==============================================================================
# Generics to replace slots in a dbFrame
# ------------------------------------------------------------------------------

setGeneric(name="bc_ids<-",      
    package="CATALYST",
    def=function(x, value) standardGeneric("bc_ids<-"))

setGeneric(name="mhl_dists<-",      
    package="CATALYST",
    def=function(x, value) standardGeneric("mhl_dists<-"))

setGeneric(name="sep_cutoffs<-", 
    package="CATALYST",
    def=function(x, value) standardGeneric("sep_cutoffs<-"))

setGeneric(name="mhl_cutoff<-",  
    package="CATALYST",
    def=function(x, value) standardGeneric("mhl_cutoff<-"))

setGeneric(name="counts<-",      
    package="CATALYST",
    def=function(x, value) standardGeneric("counts<-"))

setGeneric(name="yields<-",      
    package="CATALYST",
    def=function(x, value) standardGeneric("yields<-"))

# ==============================================================================
# Generics for debarcoding
# ------------------------------------------------------------------------------

#' @rdname assignPrelim
#' @param ... optional arguments.
#' @export
setGeneric(name="assignPrelim", 
    package="CATALYST",
    def=function(x, y, ...) standardGeneric("assignPrelim"))

#' @rdname estCutoffs
#' @param ... optional arguments.
#' @export
setGeneric(name="estCutoffs",   
    package="CATALYST",
    def=function(x, ...) standardGeneric("estCutoffs"))

#' @rdname applyCutoffs
#' @param ... optional arguments.
#' @export
setGeneric(name="applyCutoffs", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("applyCutoffs"))

#' @rdname outFrames
#' @param ... optional arguments.
#' @export
setGeneric(name="outFrames",      
    package="CATALYST", 
    def=function(x, ...) standardGeneric("outFrames"))

#' @rdname outFCS
#' @param ... optional arguments.
#' @export
setGeneric(name="outFCS",      
    package="CATALYST", 
    def=function(x, y, out_path=tempdir(), ...) standardGeneric("outFCS"))

# ==============================================================================
# Generics for plotting
# ------------------------------------------------------------------------------

#' @rdname plotYields
#' @param ... optional arguments.
#' @export
setGeneric(name="plotYields", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotYields"))

#' @rdname plotEvents
#' @param ... optional arguments.
#' @export
setGeneric(name="plotEvents", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotEvents"))

#' @rdname plotMahal
#' @param ... optional arguments.
#' @export
setGeneric(name="plotMahal", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotMahal"))

# ==============================================================================
# Generics for compensation
# ------------------------------------------------------------------------------

#' @rdname estTrim
#' @param ... optional arguments.
#' @export
setGeneric(name="estTrim", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("estTrim"))

#' @rdname computeSpillmat
#' @param ... optional arguments.
#' @export
setGeneric(name="computeSpillmat", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("computeSpillmat"))

#' @rdname plotSpillmat
#' @param ... optional arguments.
#' @export
setGeneric(name="plotSpillmat", 
    package="CATALYST",
    def=function(bc_ms, SM, ...) standardGeneric("plotSpillmat"))

#' @rdname adaptSpillmat
#' @param ... optional arguments.
#' @export
setGeneric(name="adaptSpillmat", 
    package="CATALYST",
    def=function(input_sm, out_chs) standardGeneric("adaptSpillmat"))

#' @rdname compCytof
#' @param ... optional arguments.
#' @export
setGeneric(name="compCytof", 
    package="CATALYST",
    def=function(x, y, ...) standardGeneric("compCytof"))


# ==============================================================================
# Generics to access slots in a daFrame
# ------------------------------------------------------------------------------
#' @rdname daFrame-methods
#' @export
setGeneric(name="lineage", 
    package="CATALYST", 
    def=function(x) standardGeneric("lineage"))

#' @rdname daFrame-methods
#' @export
setGeneric(name="functional",      
    package="CATALYST", 
    def=function(x) standardGeneric("functional"))

#' @rdname daFrame-methods
#' @export
setGeneric(name="sample_ids", 
    package="CATALYST", 
    def=function(x) standardGeneric("sample_ids"))

#' @rdname daFrame-methods
#' @export
setGeneric(name="conditions", 
    package="CATALYST", 
    def=function(x) standardGeneric("conditions"))

#' @rdname daFrame-methods
#' @export
setGeneric(name="cluster_codes",  
    package="CATALYST", 
    def=function(x) standardGeneric("cluster_codes"))

#' @rdname daFrame-methods
#' @export
setGeneric(name="cluster_ids",  
    package="CATALYST", 
    def=function(x) standardGeneric("cluster_ids"))

# ==============================================================================
# Generics for differential analysis
# ------------------------------------------------------------------------------

#' @rdname plotCounts
#' @param ... optional arguments.
#' @export
setGeneric(name="plotCounts", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotCounts"))

#' @rdname plotExprs
#' @param ... optional arguments.
#' @export
setGeneric(name="plotExprs", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotExprs"))

#' @rdname plotMDS
#' @param ... optional arguments.
#' @export
setGeneric(name="plotMDS", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotMDS"))

#' @rdname plotExprHeatmap
#' @param ... optional arguments.
#' @export
setGeneric(name="plotExprHeatmap", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotExprHeatmap"))

#' @rdname plotClusterHeatmap
#' @param ... optional arguments.
#' @export
setGeneric(name="plotClusterHeatmap", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotClusterHeatmap"))

#' @rdname plotMedExprs
#' @param ... optional arguments.
#' @export
setGeneric(name="plotMedExprs", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotMedExprs"))

#' @rdname plotCodes
#' @param ... optional arguments.
#' @export
setGeneric(name="plotCodes", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotCodes"))

#' @rdname plotDeltaArea
#' @param ... optional arguments.
#' @export
setGeneric(name="plotDeltaArea", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotDeltaArea"))

#' @rdname plotNRS
#' @param ... optional arguments.
#' @export
setGeneric(name="plotNRS", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotNRS"))

#' @rdname tSNE
#' @param ... optional arguments.
#' @export
setGeneric(name="tSNE", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("tSNE"))

#' @rdname plotSNE
#' @param ... optional arguments.
#' @export
setGeneric(name="plotSNE", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotSNE"))

#' @rdname mergeClusters
#' @param ... optional arguments.
#' @export
setGeneric(name="mergeClusters", 
    package="CATALYST",
    def=function(x, table, label) standardGeneric("mergeClusters"))

#' @rdname plotAbundances
#' @param ... optional arguments.
#' @export
setGeneric(name="plotAbundances", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotAbundances"))

#' @rdname diffAbundance
#' @param ... optional arguments.
#' @export
setGeneric(name="diffAbundance", 
    package="CATALYST",
    def=function(x, k, K, ...) 
        standardGeneric("diffAbundance"))

#' @rdname diffExpr
#' @param ... optional arguments.
#' @export
setGeneric(name="diffExpr", 
    package="CATALYST",
    def=function(x, k, K, ...)
        standardGeneric("diffExpr"))

#' @rdname plotDiffHeatmap
#' @param ... optional arguments.
#' @export
setGeneric(name="plotDiffHeatmap", 
    package="CATALYST",
    def=function(x, K, ...) 
        standardGeneric("plotDiffHeatmap"))