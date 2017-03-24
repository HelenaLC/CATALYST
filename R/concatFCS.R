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
#' a character string specifying the location 
#' of the FCS files to be concatinated.
#' 
#' @param y
#' can be either a \code{\link{flowSet}}, a list of \code{\link{flowFrame}}s, 
#' a character specifying the location of the FCS files to be concatinated, 
#' or a vector of FCS file names.
#' 
#' @param out_path
#' an optional character string. If specified, an FCS file 
#' of the concatinated data will be written to this location. 
#' If NULL (default), \code{concatFCS} will return a \code{\link{flowFrame}}.
#'
#' @return 
#' Returns either a \code{\link{flowFrame}} or FCS files containing 
#' the measured intensities of all provided FCS files. Files will be 
#' linked together in the order of appearance in the specified location.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom flowCore flowFrame colnames exprs 
#' @export
# ==============================================================================

setMethod(f="concatFCS",
    signature=signature(x="character"),
    definition=function(x, y=NULL) {     
        
        fcs <- list.files(x, ".fcs", full.names=TRUE, ignore.case=TRUE)
        n <- length(fcs)
        if (n < 2) {
            if (n == 0) stop("No FCS files have been found in \"", x, "\"")
            if (n == 1) stop("Only a single FCS file has been found in \"", 
                x, "\"\n  - No concatination required!")
        }
        
        ffs <- lapply(fcs, flowCore::read.FCS)
        chs <- lapply(ffs, flowCore::colnames)
        for (i in 2:n) 
            if (!isTRUE(all.equal(chs[[1]], chs[[i]])))
                stop("\n  Channel names of flowFrame ", i, 
                    " don't match those of flowFrame 1.\n ",
                    " Please make sure all FCS files share the same",
                    " measurement parameters!")
        
        n_events <- unlist(lapply(ffs, function(k) nrow(k)))
        exprs <- lapply(ffs, function(k) flowCore::exprs(k))
        time_col <- grep("time", colnames(exprs[[1]]), TRUE)
        times <- lapply(exprs, function(k) k[, time_col])
        
        new_times <- times[[1]]
        for (i in 2:n) 
            new_times <- append(new_times, times[[i]]+max(new_times))
        new_es <- do.call(rbind, exprs)
        new_es[, time_col] <- new_times
        ff <- new("flowFrame", exprs=new_es, 
            parameters=flowCore::parameters(ffs[[1]]), 
            description=flowCore::description(ffs[[1]]))
        
        if (is.null(y)) {
            return(ff)
        } else {
            suppressWarnings(flowCore::write.FCS(ff, 
                file.path(y, paste0(gsub("_\\d+.fcs", "",
                    list.files(x, ".fcs", ignore.case=TRUE)[[1]], 
                    ignore.case=TRUE), "_concat.fcs"))))
        }
    })

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
        exprs <- flowCore::fsApply(x, flowCore::exprs)
        n <- c(flowCore::fsApply(x, nrow))
        time_col <- grep("time", chs, TRUE)
        
        ts <- exprs[, time_col]
        for (i in 2:length(x)) {
            inds <- (sum(n[1:(i-1)])+1):sum(n[1:i])
            ts[inds] <- ts[inds] + max(ts[1:sum(n[1:(i-1)])])
        }
        exprs[, time_col] <- ts

        ff <- new("flowFrame", exprs=exprs, description=list())
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
