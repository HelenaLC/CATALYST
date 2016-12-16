# ==============================================================================
# Compensate CyTOF experiment
# ------------------------------------------------------------------------------

#' @rdname compCytof
#' @title Compensate CyTOF experiment
#' 
#' @description 
#' For each barcode, estimates a cutoff parameter for the 
#' distance between positive and negative barcode populations.
#'
#' @param x       
#' a \code{\link{flowFrame}} OR a character string specifying 
#' the location of FCS files that should be compensates.
#' @param y 
#' a spillover matrix.
#' @param out_path
#' a character string. If specified, compensated FCS files will be generated 
#' in this location. If \code{x} is a character string, file names will be 
#' inherited from uncompensated FCS files and given extension "_comped".
#' Defaults to NULL. 
#' 
#' @return
#' Compensates the input \code{\link{flowFrame}} OR, 
#' if \code{x} is a character string, all FCS files in the specified location.
#' 
#' @examples
#' data(ss_exp)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' spillMat <- computeSpillmat(x = re)
#' compCytof(ss_exp, spillMat)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom flowCore colnames exprs flowFrame
#' @export
# ------------------------------------------------------------------------------

setMethod(f="compCytof",
    signature=signature(x="flowFrame", y="matrix"),
    definition=function(x, y, out_path=NULL) {
        
    nms <- flowCore::colnames(x)
    ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
    ff_chs <- flowCore::colnames(x[, !is.na(ms)])
    
    # check which channels of spillover matrix are missing in flowFrame
    # and drop corresponding rows and columns
    y <- make_symetric(y)
    ex <- rownames(y)[!rownames(y) %in% ff_chs]
    if (length(ex) != 0)
        y <- y[!rownames(y) %in% ex, !colnames(y) %in% ex]
    sm_chs <- rownames(y)
    
    # check which channels of input flowFrame are not 
    # contained in spillover matrix and give warning
    add <- ff_chs[which(!ff_chs %in% sm_chs)]
    if (length(add) != 0) {
        new_mets <- gsub("[[:digit:]]+Di", "", add)
        old_ms <- as.numeric(regmatches(sm_chs, gregexpr("[0-9]+", sm_chs)))
        new_ms <- as.numeric(regmatches(add, gregexpr("[0-9]+", add)))
        nms <- c(sm_chs, add)[order(c(old_ms, new_ms))]
        spill_cols <- get_spill_cols(new_ms, new_mets)
        
        message("WARNING: Compensation is likely to be inaccurate.\n",
            "         Spill values for the following interactions\n",
            "         have not been estimated and are missing:")
        for (i in seq_along(add)) {
            if (length(spill_cols[[i]]) == 0) next
            cat(add[i], "->", paste(nms[spill_cols[[i]]], collapse=", "), "\n")
        }
        # add them into the matrix
        sm <- diag(length(nms))
        rownames(sm) <- colnames(sm) <- nms
        sm[sm_chs, sm_chs] <- y
    }
    compensate(x, y)
    })

setMethod(f="compCytof",
    signature=signature(x="character", y="matrix"),
    definition=function(x, y, out_path=NULL) {
        if (!file.exists(x))
            stop("x is nor a flowFrame nor a valid file/folder path.")
        fcs <- list.files(path=x, patter=".fcs", full.names=TRUE)
        if (length(fcs) == 0)
            stop("No FCS files found in specified location.")
        ffs <- lapply(fcs, flowCore::read.FCS)
        out_nms <- paste0(gsub(".fcs", "", fcs), "_comped.fcs")
        for (i in seq_along(ffs))
            write.FCS(compCytof(ffs[[i]], y), out_nms[i])
        })






