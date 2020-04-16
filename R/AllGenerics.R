#' @export
setGeneric("ei", 
    function(x) standardGeneric("ei"))

#' @export
setGeneric("n_cells", 
    function(x) standardGeneric("n_cells"))

#' @export
setGeneric("channels", 
    function(x) standardGeneric("channels"))

#' @export
setGeneric("marker_classes", 
    function(x) standardGeneric("marker_classes"))

#' @export
setGeneric("type_markers", 
    function(x) standardGeneric("type_markers"))

#' @export
setGeneric("state_markers", 
    function(x) standardGeneric("state_markers"))

#' @export
setGeneric("sample_ids", 
    function(x) standardGeneric("sample_ids"))

#' @export
setGeneric("cluster_ids",  
    function(x, k) 
        standardGeneric("cluster_ids"))

#' @export
setGeneric("cluster_codes",  
    function(x) standardGeneric("cluster_codes"))

#' @export
setGeneric("delta_area",  
    function(x) standardGeneric("delta_area"))