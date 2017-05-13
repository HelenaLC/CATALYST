# ==============================================================================
# flowFrames from dbFrame
# ------------------------------------------------------------------------------
#' @rdname outFrames
#' @title Population-wise \code{flowFrame}s from a \code{dbFrame}
#' @description Returns a \code{flowSet} or list of \code{flowFrame}s from a 
#' \code{\link{dbFrame}}. Each \code{flowFrame} will contain the subset of 
#' events that have been assigned to the same ID. 
#'
#' @param x 
#' a \code{\link{dbFrame}}.
#' @param return
#' \code{"flowSet"} or \code{"list"}. Specifies the output type.
#' @param which
#' Specifies which barcode(s) to include. \code{"assigned"} (if the population
#' of unassigned events should be excluded), \code{"all"} (if the latter should
#' be included), or a numeric or character specifying a subset of populations. 
#' Valid values are IDs that occur as row names in the \code{bc_key} of the 
#' supplied \code{\link{dbFrame}}. Defaults to \code{"assigned"}.
#' 
#' @details 
#' Creates a separate \code{\link{flowFrame}} for each barcode population 
#' and, if desired, the population of unassigned events.
#' 
#' @return a \code{\link{flowSet}} or list of \code{\link{flowFrame}}s.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' outFfs(x = re, return = "list", which = c("B1", "D4"))
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom flowCore flowFrame flowSet
#' @export
# ------------------------------------------------------------------------------

setMethod(f="outFrames",       
    signature=signature(x="dbFrame"), 
    definition=function(x, return="flowSet", which="assigned") {
        
        if (!return %in% c("flowSet", "list"))
            stop("Invalid 'return' argument. Valid values are 
                \"flowSet\" or \"list\" for a list of flowFrames.")
        
        ids <- NULL
        if (length(which) == 1) {
            if (which == "assigned") {
                ids <- sort(unique(bc_ids(x)))
                ids <- ids[ids != 0]
            } else if (which == "all") {
                ids <- sort(unique(bc_ids(x)))
            }
        }
        if (is.null(ids)) {
            test <- which %in% c(0, rownames(bc_key(x)))
            if (sum(test) != length(which))
                stop("Valid values for 'which' are IDs that occur as row names
                    in the\n", " 'bc_key' slot of the supplied 'dbFrame', or 0 
                    for unassigned events.")
            ids <- which
        }
        
        n <- length(ids)
        ffs <- vector("list", n)
        
        for (i in seq_along(ids)) {
            inds <- bc_ids(x) == ids[i]
            if (sum(inds) < 2) 
                next
            ffs[[i]] <- new("flowFrame", exprs=exprs(x)[inds, ])
            flowCore::description(ffs[[i]])$GUID <- ids[i]
        }
        
        empty <- which(sapply(ffs, is.null))
        nSkipped <- length(empty)
        if (nSkipped > 0) {
            ffs <- ffs[-empty] 
            if (nSkipped == n) {
                stop("No or less than 2 events have been assigned 
                    to the specified population(s).")
            } else if (nSkipped == 1) {
                warning("ID ", ids[empty], 
                    " has less than 2 event assignments 
                    and has been skipped.")
            } else {
                warning("IDs ", paste(ids[empty], collapse=", "), 
                    " have less than 2 event assignments 
                    and have been skipped.")
            }
            
            }
        
        if (return == "flowSet") {
            set <- as(ffs, "flowSet")
            flowCore::sampleNames(set) <- sapply(ffs, 
                function(i) description(i)$GUID)
            set
        } else {
            ffs
        }
    })