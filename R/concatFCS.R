# ==============================================================================
# concatinates all FCS files in the specified location
# ------------------------------------------------------------------------------

#' @rdname concatFCS
#' @title FCS file concatination
#' 
#' @description 
#' Concatinates all input data.
#'
#' @param x 
#' can be either a \code{\link{flowSet}}, a list of \code{\link{flowFrame}}s, 
#' a character specifying the location of the FCS files to be concatinated, 
#' or a vector of FCS file names.
#' 
#' @param out_path
#' an optional character string. If specified, an FCS file 
#' of the concatinated data will be written to this location. 
#' If NULL (default), a \code{\link{flowFrame}} will be returned.
#'
#' @return 
#' a \code{\link{flowFrame}} containing measurement intensities 
#' of all input data or a character of the FCS file name.
#' 
#' @examples
#' data(raw_data)
#' concatFCS(raw_data)
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom flowCore colnames exprs flowFrame flowSet fsApply isFCSfile
#' @importFrom methods as
#' @export

# ------------------------------------------------------------------------------

setMethod(f="concatFCS",
    signature=signature(x="flowSet"),
    definition=function(x, out_path=NULL) {
        
        chs <- flowCore::colnames(x)
        es <- flowCore::fsApply(x, flowCore::exprs)
        n <- flowCore::fsApply(x, nrow)
        t <- grep("time", chs, TRUE)

        start <- c(1, cumsum(n)+1)
        end <- start[-1]-1
        for (i in seq_along(x)[-1]) {
            inds <- start[i]:end[i]
            es[inds, t] <- es[inds, t] + es[end[i-1], t]
        }
        
        ff <- new("flowFrame", exprs=es, description=list())
        flowCore::parameters(ff)$desc <- flowCore::parameters(x[[1]])$desc
        
        nm <- gsub("_\\d+.fcs", "_concat.fcs", 
            flowCore::description(x[[1]])$ORIGINALGUID, TRUE)
        flowCore::description(ff)[c("GUID", "ORIGINALGUID")] <- 
            flowCore::identifier(ff) <- nm
        
        if (is.null(out_path)) 
            return(ff)
        suppressWarnings(flowCore::write.FCS(ff, file.path(out_path, nm)))
    })

# ------------------------------------------------------------------------------

#' @rdname concatFCS
setMethod(f="concatFCS",
    signature=signature(x="character"),
    definition=function(x, out_path=NULL) {
        if (length(x) == 1) {
            fcs <- list.files(x, ".fcs", full.names=TRUE, ignore.case=TRUE)
        } else {
            fcs <- x
        } 
        n <- length(fcs)
        if (sum(flowCore::isFCSfile(fcs)) != n) 
            stop("Not all files in ", x, " are valid FCS files.")   
        if (n < 2) {
            if (n == 0) 
                stop("No FCS files have been found in \"", x, "\"")
            if (n == 1) 
                stop("Only a single FCS file has been found in \"", x,"\"") 
        }
        concatFCS(lapply(fcs, flowCore::read.FCS), out_path)
    })

# ------------------------------------------------------------------------------

#' @rdname concatFCS
setMethod(f="concatFCS",
    signature=signature(x="list"),
    definition=function(x, out_path=NULL) {
        if (length(x) == 1) 
            stop("Only a single flowFrame has been provided.")
        concatFCS(as(x, "flowSet"), out_path)
    })
