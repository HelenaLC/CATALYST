# ==============================================================================
# Generics to access slots in a dbFrame
# ------------------------------------------------------------------------------

#' @rdname dbFrame-methods
#' @export
setGeneric(name="bc_key", 
    package="CATALYST", 
    def=function(object) standardGeneric("bc_key"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="bc_ids",      
    package="CATALYST", 
    def=function(object) standardGeneric("bc_ids"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="deltas",      
    package="CATALYST", 
    def=function(object) standardGeneric("deltas"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="normed_bcs",  
    package="CATALYST", 
    def=function(object) standardGeneric("normed_bcs"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="sep_cutoffs", 
    package="CATALYST", 
    def=function(object) standardGeneric("sep_cutoffs"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="mhl_cutoff",  
    package="CATALYST", 
    def=function(object) standardGeneric("mhl_cutoff"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="counts",      
    package="CATALYST", 
    def=function(object) standardGeneric("counts"))

#' @rdname dbFrame-methods
#' @export
setGeneric(name="yields",      
    package="CATALYST", 
    def=function(object) standardGeneric("yields"))

# ==============================================================================
# Generics to replace slots in a dbFrame
# ------------------------------------------------------------------------------

# required by estCutoffs()
setGeneric(name="sep_cutoffs<-", 
    package="CATALYST",
    def=function(object, value) standardGeneric("sep_cutoffs<-"))

setGeneric(name="counts<-",      
    package="CATALYST",
    def=function(object, value) standardGeneric("counts<-"))

setGeneric(name="yields<-",      
    package="CATALYST",
    def=function(object, value) standardGeneric("yields<-"))

# required by applyCutoffs()
setGeneric(name="bc_ids<-",      
    package="CATALYST",
    def=function(object, value) standardGeneric("bc_ids<-"))

setGeneric(name="mhl_cutoff<-",  
    package="CATALYST",
    def=function(object, value) standardGeneric("mhl_cutoff<-"))

# ==============================================================================
# Generics for debarcoding
# ------------------------------------------------------------------------------
#' @rdname assignPrelim
#' @param ... further optional arguments.
#' @export
setGeneric(name="assignPrelim", 
    package="CATALYST",
    def=function(x, y, ...) 
        standardGeneric("assignPrelim"))

#' @rdname estCutoffs
#' @param ... further optional arguments.
#' @export
setGeneric(name="estCutoffs",   
    package="CATALYST",
    def=function(x, ...) 
        standardGeneric("estCutoffs"))

#' @rdname applyCutoffs
#' @param ... further optional arguments.
#' @export
setGeneric(name="applyCutoffs", 
    package="CATALYST",
    def=function(x, ...) 
        standardGeneric("applyCutoffs"))

#' @rdname outFCS
#' @param ... further optional arguments.
#' @export
setGeneric(name="outFCS",      
    package="CATALYST", 
    def=function(x, out_path) standardGeneric("outFCS"))

#' @rdname debarcode
#' @param ... further optional arguments.
#' @export
setGeneric(name="debarcode",
    package="CATALYST",
    def=function(x, y, out_path, ...) 
        standardGeneric("debarcode"))

# ==============================================================================
# Generics for compensation
# ------------------------------------------------------------------------------
#' @rdname estTrim
#' @export
setGeneric(name="estTrim", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("estTrim"))

#' @rdname computeCompmat
#' @param ... further optional arguments.
#' @export
setGeneric(name="computeCompmat", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("computeCompmat"))

# ==============================================================================
# Generics for plotting
# ------------------------------------------------------------------------------
#' @rdname plotYields
#' @param ... further optional arguments.
#' @export
setGeneric(name="plotYields", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotYields"))

#' @rdname plotEvents
#' @param ... further optional arguments.
#' @export
setGeneric(name="plotEvents", 
    package="CATALYST",
    def=function(x, ...) standardGeneric("plotEvents"))

