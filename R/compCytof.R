# ==============================================================================
# Compensate CyTOF experiment
# ------------------------------------------------------------------------------

#' @rdname compCytof
#' @title Compensate CyTOF experiment
#' 
#' @description 
#' Compensates a mass spectrometry based experiment using a provided spillover
#' matrix, assuming a linear spillover in the experiment.
#' 
#' If the spillover matrix does not contain all the same columns than the experiment,
#' it will be adapted according to the following rules:
#' - non metal columns present in the experiment but not in the psillover matrix will
#'   be added such that they do neiter receive nor emit spillover
#'   -> Exception: if the added metal has a mass equal to amass already present in the
#'      spillover matrix, it will receive (but not emit) spillover according to the present
#'      metal with the same mass.
#'   -> If an added metal could potentially receive spillover, as it is the same metal type, has a
#'      mass of +- 1 or +16 of a metal present in the experiment a warning will be issued,
#'      as there could be a possible spillover interaction missed in the experiment, leading
#'      potentially to faulty compensation.
#'      
#' - columns present in the spillover matrix but not in the experiment will be removed from
#'   the spillover matrix.
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
#' @details
#' Compensates the input \code{\link{flowFrame}} OR, 
#' if \code{x} is a character string, all FCS files in the specified location.
#' 
#' @return 
#' If \code{out_path} is NULL (default), returns a \code{\link{flowFrame}} 
#' containing the compensated measurement data. Else, compensated data 
#' will be  written to the specified location as FCS 3.0 standard files. 
#' 
#' @examples
#' data(ss_beads)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_beads, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' spillMat <- computeSpillmat(x = re)
#' compCytof(ss_beads, spillMat)
#'
#' @author 
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' and Vito Zanotelli \email{vito.zanotelli@uzh.ch}
#' @importFrom flowCore flowFrame colnames exprs compensate
#' @export
# ------------------------------------------------------------------------------

setMethod(f="compCytof",
    signature=signature(x="flowFrame", y="matrix"),
    definition=function(x, y, out_path=NULL) {
        
        # check validity of input spillover matrix
        if (any(y < 0))
            stop("\nThe supplied spillover matrix is invalid ",
                "as it contains negative entries.\n",
                "Valid spillvalues are non-negative and mustn't exceed 1.")
        if (any(y > 1))
            stop("\nThe supplied spillover matrix is invalid ",
                "as it contains entries greater than 1.\n",
                "Valid spillvalues are non-negative and mustn't exceed 1.")
        
        nms <- flowCore::colnames(x)
        ms <- gsub("[[:alpha:][:punct:]]", "", nms)
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
                idx <- which(ms == new_ms[i] & all_mets == new_mets[i])
                if ( length(idx) > 0) {
                    if (first) {
                        message("WARNING: ",
                            "Compensation is likely to be inaccurate.\n",
                            "         ",
                            "Spill values for the following interactions\n",
                            "         ",
                            "have not been estimated:")
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
                # add the spillover
                sm[rownames(y), sm_col] <- y[,y_col[sm_col_ms]]
                for (m in unique(sm_col_ms)){
                    mfil = ms == m
                    # set the spillover between channels of the same mass to 0
                    # other wise the linear system can get singular.
                    # the diagonal elements will be set to 1 lateron again
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

#' @rdname compCytof
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


