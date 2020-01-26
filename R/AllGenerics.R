# ==============================================================================
# Generics to access slots in a SCE constructed with prepData()
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