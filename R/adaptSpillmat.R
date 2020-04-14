#' @rdname adaptSpillmat
#' @title Adapt spillover matrix
#' 
#' @description 
#' This helper function adapts the columns of a provided spillover matrix 
#' such that it is compatible with data having the column names provided.
#'
#' @param x a previously calculated spillover matrix.
#' @param out_chs the column names that the prepared output 
#'   spillover matrix should have. Numeric names as well as names 
#'   of the form MetalMass(Di), e.g. Ir191, Ir191Di or Ir191(Di), 
#'   will be interpreted as masses with associated metals.
#' @param isotope_list named list. Used to validate the input spillover matrix.
#'   Names should be metals; list elements numeric vectors of their isotopes.
#'   See \code{\link{isotope_list}} for the list of isotopes used by default.
#' @param verbose logical. Should warnings about possibly 
#'   inaccurate spillover estimates be printed to the console?
#'  
#' @details The rules how the spillover matrix is adapted 
#' are explained in \code{\link{compCytof}}. 
#' 
#' @return An adapted spillover matrix with 
#' column and row names according to \code{out_chs}.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch} 
#' & Vito RT Zanotelli
#' 
#' @examples
#' # estimate spillover matrix from 
#' # single-stained control samples
#' data(ss_exp)
#' sce <- prepData(ss_exp)
#' bc_ms <- c(139, 141:156, 158:176)
#' sce <- assignPrelim(sce, bc_ms, verbose = FALSE)
#' sce <- applyCutoffs(estCutoffs(sce))
#' sce <- computeSpillmat(sce)
#' 
#' library(SingleCellExperiment)
#' sm1 <- metadata(sce)$spillover_matrix
#' sm2 <- adaptSpillmat(sm1, rownames(sce), verbose = FALSE)
#' all(dim(sm2) == ncol(sm1))
#' 
#' @export

adaptSpillmat <- function(x, out_chs, 
    isotope_list = CATALYST::isotope_list, 
    verbose = TRUE) {
    # check validity of input spillover matrix
    sm <- .check_sm(x, isotope_list)
    # get the output names, metals and masses
    n <- length(out_chs)
    out_masses <- .get_ms_from_chs(out_chs)
    out_metalchs <- out_chs[!is.na(out_masses)]
    input_sm_chs_col <- colnames(x)
    input_sm_chs_row <- rownames(x)
    input_sm_chs <- unique(c(input_sm_chs_row, input_sm_chs_col))
    
    # copy the existing spillover information into a new spillover matrix
    sm <- matrix(diag(n), n, n, dimnames=list(out_chs, out_chs))
    sm_preexisting_col <- input_sm_chs_col[input_sm_chs_col %in% out_chs]
    sm_preexisting_row <- input_sm_chs_row[input_sm_chs_row %in% out_chs]
    sm[sm_preexisting_row, sm_preexisting_col] <- 
        x[sm_preexisting_row, sm_preexisting_col]
    
    if (verbose) .warn_new_intearctions(out_metalchs, x)
    
    # check for new channels
    new_metalchs <- out_metalchs[!out_metalchs %in% input_sm_chs]
    new_masses <- .get_ms_from_chs(new_metalchs)
    
    old_receiving_masses <- .get_ms_from_chs(input_sm_chs_col)
    
    test <- (length(new_metalchs) != 0) && 
        (any(inds <- old_receiving_masses %in% new_masses))
    if (test) {
        # check if any new masses were already present in the old masses
        # and add them to receive spillover according to the old masses
        
        # get the channels that correspond to the old_masses 
        # that have an aditional metal with the same weight
        y_col <- input_sm_chs_col[inds]
        names(y_col) <- as.character(old_receiving_masses[inds])
        # get all columns that are part of the affected masses
        fil <- out_masses %in% old_receiving_masses[inds]
        sm_col <- out_chs[fil]
        sm_col_ms <- as.character(out_masses[fil])
        # add the spillover
        old_rowchs <- out_chs[out_chs %in% input_sm_chs_row]
        sm[old_rowchs, sm_col] <- x[old_rowchs, y_col[sm_col_ms]]
        for (m in unique(sm_col_ms)) {
            mfil <- out_masses == m
            # set the spillover between channels of the same mass to 0
            # otherwise the linear system can get singular
            # diagonal elements will be set to 1 again later on
            sm[mfil, mfil] <- 0
        }
    }
    # assure diagonal is all 1
    diag(sm) <- 1
    return(sm)  
}