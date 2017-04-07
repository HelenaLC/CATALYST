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
#' If the spillover matrix (SM) does not contain the same set of columns as 
#' the input experiment, it will be adapted according to the following rules:
#' \enumerate{
#' \item{columns present in the SM but not in the input data 
#' will be removed from it}
#' \item{non-metal columns present in the input but not in the SM 
#' will be added such that they do neither receive nor cause spill}
#' \item{metal columns that have the same mass as a channel present in the SM 
#' will receive (but not emit) spillover according to that channel}
#' \item{if an added channel could potentially receive spillover (as it has 
#' +/-1M or +16M of, or is of the same metal type as another channel measured), 
#' a warning will be issued as there could be spillover interactions that have
#' been missed and may lead to faulty compensation}
#' }
#' 
#' @return 
#' Compensates the input \code{\link{flowFrame}} or, 
#' if \code{x} is a character string, all FCS files in the specified location. 
#' If \code{out_path=NULL} (the default), returns a \code{\link{flowFrame}} 
#' containing the compensated data. Otherwise, compensated data will be written 
#' to the specified location as FCS 3.0 standard files. 
#' 
#' @examples
#' # get single-stained control samples
#' # get single-stained control samples
#' data(ss_exp)
#' 
#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#' 
#' # debarcode
#' re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' spillMat <- computeSpillmat(x = re)
#' compCytof(x = ss_exp, y = spillMat)
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
        # if (any(y < 0))
        #     stop("\nThe supplied spillover matrix is invalid ",
        #         "as it contains negative entries.\n",
        #         "Valid spillvalues are non-negative and mustn't exceed 1.")
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
            all_mets <- gsub("[[:digit:]]+Di", "", nms)
            spill_cols <- get_spill_cols(ms, all_mets)
            
            first = TRUE
            for (i in seq_along(new_ms)) {
                idx <- which(ms == new_ms[i] & all_mets == new_mets[i])
                if (length(idx) > 0) {
                    if (first) {
                        message("WARNING: ",
                            "Compensation is likely to be inaccurate.\n",
                            "         ",
                            "Spill values for the following interactions\n",
                            "         ",
                            "have not been estimated:")
                        first=FALSE
                    }
                    message(nms[idx], " -> ", paste(
                        nms[spill_cols[[idx]]], collapse=", "))
                }
            }
        }
        
        # add them into the matrix
        sm <- diag(length(nms))
        rownames(sm) <- colnames(sm) <- nms
        sl_sm_cols <- sm_cols[sm_cols %in% ff_chs]
        sm[sm_chs, sl_sm_cols] <- y[sm_chs, sl_sm_cols]
        
        test <- (length(add) != 0) && (any(inds <- old_ms %in% new_ms))
        if (test) {
            # check if any new masses were already present in the old masses
            # and add them to receive spillover according to the old masses
            
            # get the channels that correspond to the old_masses 
            # that have an aditional metal with the same weight
            y_col <- sm_chs[inds]
            names(y_col) <- as.character(old_ms[inds])
            # get all columns that are part of the affected masses
            fil <- ms %in% old_ms[inds]
            sm_col <- nms[fil]
            sm_col_ms <- as.character(ms[fil])
            # add the spillover
            sm[rownames(y), sm_col] <- y[, y_col[sm_col_ms]]
            for (m in unique(sm_col_ms)){
                mfil <- ms == m
                # set the spillover between channels of the same mass to 0
                # otherwise the linear system can get singular
                # diagonal elements will be set to 1 again later on
                sm[mfil, mfil] <- 0
            }
        }
        
        # check which channels of spillover matrix are missing in flowFrame
        # and drop corresponding rows and columns
        ex <- rownames(sm)[!rownames(sm) %in% ff_chs]
        if (length(ex) != 0)
            sm <- sm[!rownames(sm) %in% ex, !colnames(sm) %in% ex]
        
        # assure diagonal is all 1
        diag(sm) <- 1
        
        comped <- flowCore::compensate(x, sm)
        if (!is.null(out_path)) {
            nm <- flowCore::identifier(x)
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
            stop("x is neither a flowFrame nor a valid file/folder path.")
        fcs <- list.files(x, ".fcs", full.names=TRUE)
        if (length(fcs) == 0)
            stop("No FCS files found in specified location.")
        ffs <- lapply(fcs, flowCore::read.FCS)
        
        if (is.null(out_path)) {
            out_nms <- gsub(".fcs$", "_comped.fcs", fcs)
            lapply(ffs, function(i) compCytof(i, y))
        } else {
            out_nms <- file.path(out_path, gsub(".fcs", 
                "_comped.fcs", list.files(x, ".fcs")))
            for (i in seq_along(ffs))
                suppressWarnings(flowCore::write.FCS(
                    compCytof(ffs[[i]], y), out_nms[i]))
        }
    })
