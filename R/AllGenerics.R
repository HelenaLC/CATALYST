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

#' @rdname SCE-utils
#' @export
setGeneric("ei", 
    function(x) standardGeneric("ei"))

#' @rdname SCE-utils
#' @export
setGeneric("n_cells", 
    function(x) standardGeneric("n_cells"))

#' @rdname SCE-utils
#' @export
setGeneric("marker_classes", 
    function(x) standardGeneric("marker_classes"))

#' @rdname SCE-utils
#' @export
setGeneric("type_markers", 
    function(x) standardGeneric("type_markers"))

#' @rdname SCE-utils
#' @export
setGeneric("state_markers",      
    function(x) standardGeneric("state_markers"))

#' @rdname SCE-utils
#' @export
setGeneric("sample_ids", 
    function(x) standardGeneric("sample_ids"))

#' @rdname SCE-utils
#' @export
setGeneric("cluster_ids",  
    function(x, k) standardGeneric("cluster_ids"))

#' @rdname SCE-utils
#' @export
setGeneric("cluster_codes",  
    function(x) standardGeneric("cluster_codes"))

# ==============================================================================
# Generics for differential analysis
# ------------------------------------------------------------------------------

#' @rdname guessPanel
#' @param ... optional arguments.
#' @export
setGeneric("guessPanel", 
    function(x, ...) standardGeneric("guessPanel"))