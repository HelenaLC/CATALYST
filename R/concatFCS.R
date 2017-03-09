# ==============================================================================
# concatinates all FCS files in the specified location
# ------------------------------------------------------------------------------

#' @rdname conactFCS
#' @title FCS file concatination
#' 
#' @description 
#' Concatinates all FCS files in the specified location
#'
#' @param x 
#' a character string specifying the location 
#' of the FCS files to be concatinated.
#' @param y
#' an optional character string. If specified, an FCS file 
#' of the concatinated data will be written to this location. 
#' If NULL (default), \code{concatFCS} will return a \code{\link{flowFrame}}.
#'
#' @return 
#' Returns either a \code{\link{flowFrame}} or FCS files containing 
#' the measured intensities of all provided FCS files. Files will be 
#' linked together in the order of appearance in the specified location.
#' 
#' @examples
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom flowCore flowFrame colnames exprs 
# ==============================================================================
 
setMethod(f="concatFCS",
    signature=signature(x="character"),
    definition=function(x, y=NULL) {     
    
    fcs <- list.files(x, ".fcs", full.names=TRUE, ignore.case=TRUE)
    n <- length(fcs)
    if (n == 0) stop("No FCS files have been found in \"", x, "\"")
    if (n == 1) stop("Only a single FCS file has been found in \"", 
        x, "\"\n  - No concatination required!")
    
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
        parameters=parameters(ffs[[1]]), 
        description=description(ffs[[1]]))
    if (is.null(y)) {
        return(ff)
    } else {
        suppressWarnings(write.FCS(ff, file.path(y, paste0(gsub(".fcs", "", 
            list.files(x, ".fcs", ignore.case=TRUE)[[1]], ignore.case=TRUE), 
            "_concat.fcs"))))
    }
})