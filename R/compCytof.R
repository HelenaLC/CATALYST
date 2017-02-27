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
#' @importFrom flowCore flowFrame colnames exprs compensate
#' @export
# ------------------------------------------------------------------------------

setMethod(f="compCytof",
    signature=signature(x="flowFrame", y="matrix"),
    definition=function(x, y, out_path=NULL) {
        
    nms <- flowCore::colnames(x)
    ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
    ff_chs <- flowCore::colnames(x[, !is.na(ms)])
    sm_chs <- rownames(y)
    sm_cols <- colnames(y)
    y <- make_symetric(y)
    
    # check which channels of input flowFrame are not 
    # contained in spillover matrix and give warning
    add <- ff_chs[(!ff_chs %in% sm_chs)]
    if (length(add) != 0) {
        new_mets <- gsub("[[:digit:]]+Di", "", add)
        old_ms <- as.numeric(regmatches(sm_chs, gregexpr("[0-9]+", sm_chs)))
        new_ms <- as.numeric(regmatches(add, gregexpr("[0-9]+", add)))
        ms <- c(old_ms, new_ms)
        o <- order(ms)
        ms <- ms[o]
        nms <- c(sm_chs, add)[o]
        # get the potential spillover interactions 
        all_mets = gsub("[[:digit:]]+Di", "", nms)
        spill_cols <- get_spill_cols(ms, all_mets)
        
        first = TRUE
        for (i in seq_along(new_ms)) {
            idx = which(ms == new_ms[i] & all_mets == new_mets[i])
            if ( length(spill_cols[[idx]]) > 0) {
                if (first) {
                    message("WARNING: Compensation is likely to be inaccurate.\n",
                            "         Spill values for the following interactions\n",
                            "         have not been estimated:")
                    first = FALSE
                }
                cat(nms[idx], "->", 
                    paste(nms[spill_cols[[idx]]], collapse=", "), "\n")
            }
        }
    }
    
    # add them into the matrix
    sm <- diag(length(nms))
    rownames(sm) <- colnames(sm) <- nms
    sl_sm_cols = sm_cols[sm_cols %in% ff_chs]
    sm[sm_chs, sl_sm_cols] <- y[sm_chs, sl_sm_cols]
    
    if (length(add) != 0) {
        if (any(ind <- old_ms %in% new_ms)) {
            # check if any new masses were already present in the old masses
            # and add them to receive spillover according to the old masses
            
            # get the channels that correspond to the old_masses 
            # that have an aditional metal with the same weight
            y_col <- sm_chs[ind]
            names(y_col) <- sapply(old_ms[ind], as.character)
            # get all columns that are part of the affected masses
            fil = ms %in% old_ms[ind]
            sm_col <- nms[fil]
            sm_col_ms <-sapply(ms[fil], as.character)
            sm[rownames(y), sm_col] <- y[,y_col[sm_col_ms]]
            for (m in unique(sm_col_ms)){
                mfil = ms == m
                sm[mfil,mfil] <- 0
            }
        }
    }

    # check which channels of spillover matrix are missing in flowFrame
    # and drop corresponding rows and columns
    ex <- rownames(sm)[!rownames(sm) %in% ff_chs]
    if (length(ex) != 0)
        sm <- sm[!rownames(sm) %in% ex, !colnames(sm) %in% ex]
    
    # make sure the diagonal is all 1
    diag(sm) <- 1
    
    comped <- flowCore::compensate(x, sm)
    if (!is.null(out_path)) {
        nm <- deparse(substitute(ss_exp))
        suppressWarnings(flowCore::write.FCS(comped,
            file.path(out_path, paste0(nm, "_comped.fcs"))))
    } else {
        comped
    }
    })

# ------------------------------------------------------------------------------

setMethod(f="compCytof",
    signature=signature(x="character", y="matrix"),
    definition=function(x, y, out_path=NULL) {
        if (!file.exists(x))
            stop("x is not a flowFrame nor a valid file/folder path.")
        fcs <- list.files(path=x, pattern=".fcs", full.names=TRUE)
        if (length(fcs) == 0)
            stop("No FCS files found in specified location.")
        ffs <- lapply(fcs, flowCore::read.FCS)
        if (is.null(out_path)) {
            out_nms <- paste0(gsub(".fcs", "", fcs), "_comped.fcs")
            for (i in seq_along(ffs))
                suppressWarnings(flowCore::write.FCS(
                    compCytof(ffs[[i]], y), out_nms[i]))
        } else {
            out_nms <- file.path(out_path, list.files(path=x, pattern=".fcs"))
            for (i in seq_along(ffs))
                suppressWarnings(flowCore::write.FCS(
                    compCytof(ffs[[i]], y), out_nms[i]))
        }
        })


get_warnings <- function() {
    print("Warning") 
    return("Data")
}



