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

#' @rdname outFCS
#' @param ... optional arguments.
#' @export
setGeneric(name="outFCS",      
    package="CATALYST", 
    def=function(x, out_path=tempdir(), ...) standardGeneric("outFCS"),
    signature="x")

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

#' @rdname compCytof
#' @param ... optional arguments.
#' @export
setGeneric(name="compCytof", 
    package="CATALYST",
    def=function(x, y, ...) standardGeneric("compCytof"))