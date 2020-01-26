# ==============================================================================
# get spillover columns
# ------------------------------------------------------------------------------
.get_spill_chs <- function(ms, mets, l=CATALYST::isotope_list) {
    ms <- as.numeric(ms)
    lapply(seq_along(ms), function(i) {
        p1 <- m1 <- ox <- iso <- NULL
        if ((ms[i] + 1)  %in% ms) p1 <- which(ms == (ms[i] + 1))
        if ((ms[i] - 1)  %in% ms) m1 <- which(ms == (ms[i] - 1)) 
        if ((ms[i] + 16) %in% ms) ox <- which(ms == (ms[i] + 16))
        iso <- l[[mets[i]]]
        iso <- which(ms %in% iso[iso != ms[i]])
        unique(c(m1, p1, iso, ox))
    })
}

# ==============================================================================
# compute channel i to j spill
# ------------------------------------------------------------------------------
#' @importFrom stats median
.get_sij <- function(pos_i, neg_i, pos_j, neg_j, method, trim) {
    if (length(neg_i) == 0) neg_i <- 0
    if (length(neg_j) == 0) neg_j <- 0
    if (method == "default") {
        bg_i <- mean(neg_i, trim = 0.1)
        bg_j <- mean(neg_j, trim = 0.1)
        spiller  <- pos_i - bg_i
        receiver <- pos_j - bg_j
    } else if (method == "classic") {
        spiller  <- mean(pos_i, trim) - mean(neg_i, trim)
        receiver <- mean(pos_j, trim) - mean(neg_j, trim)
    } else {
        stop("'method = ", method, "' is not a valid option.")
    }
    receiver[receiver < 0] <- 0
    spiller [spiller  < 0] <- 0
    sij <- receiver / spiller
    sij[is.infinite(sij) | is.na(sij)] <- 0
    median(sij)
}

# ==============================================================================
# make spillover matrix symmetrical
# ------------------------------------------------------------------------------
.make_symetric <- function(x) {
    x <- as.matrix(x)
    dims <- dim(x)
    i <- which.max(dim(x))
    nms <- dimnames(x)[[i]]
    M <- matrix(0, dims[i], dims[i])
    rownames(M) <- colnames(M) <- nms
    M[rownames(x), colnames(x)] <- as.matrix(x)
    M
}

# ==============================================================================
# check validity of input spillover matrix in compCytof()
# ------------------------------------------------------------------------------
.check_sm <- function(sm, l=CATALYST::isotope_list) {
    if (any(sm < 0))
        stop("\nThe supplied spillover matrix is invalid ",
            "as it contains negative entries.\n",
            "Valid spill values are non-negative and mustn't exceed 1.")
    if (any(sm > 1))
        stop("\nThe supplied spillover matrix is invalid ",
            "as it contains entries greater than 1.\n",
            "Valid spill values are non-negative and mustn't exceed 1.")
    chs <- colnames(sm)[which(colnames(sm) %in% rownames(sm))]
    sii <- sm[cbind(chs, chs)]
    if (any(sii != 1))
        stop("\nThe supplied spillover matrix is invalid ",
            "as its diagonal contains entries != 1.\n")
    if (!all(rownames(sm) %in% colnames(sm)))
        stop("\nThe supplied spillover matrix seems to be invalid.\n",
            "All spill channels must appear as receiving channels:\n",
            "'all(rownames(sm) %in% colnames(sm))' should return TRUE.")
    isos <-  paste0(gsub("[0-9]", "", names(unlist(l))), as.numeric(unlist(l)))
    valid <- vapply(dimnames(sm), function(chs) {
        ms <- .get_ms_from_chs(chs)
        mets <- .get_mets_from_chs(chs)
        all(paste0(mets, ms) %in% isos) 
    }, logical(1))
    if (any(!valid)) 
        stop("\nThe supplied spillover matrix seems to be invalid.\n",
            "All isotopes should appear in `", deparse(substitute(l)), "`.")
    sm[, colSums(sm) != 0]
}

# ==============================================================================
# Helper functions to get mass and metal from a channel name
# ------------------------------------------------------------------------------
.get_ms_from_chs <- function(chs)
    as.numeric(gsub("[[:punct:][:alpha:]]", "", chs))

.get_mets_from_chs <- function(chs)
    gsub("([[:punct:]]*)([[:digit:]]*)((Di)|(Dd))*", "", chs)

# ==============================================================================
# This function compares a list of spillover channels to
# a provided spillover matrix and provides warnings for cases
# where the spillovermatrix has no information about a potential
# expected spillover among the new channels.
# ------------------------------------------------------------------------------
.warn_new_intearctions <- function(chs_new, sm) {
    chs <- c(list(chs_new), dimnames(sm))
    names(chs) <- c("new", "emitting", "receiving")
    
    # get masses & metals from channel names
    ms <- lapply(chs, .get_ms_from_chs)
    mets <- lapply(chs, .get_mets_from_chs)
    
    # get the potential mass channels a channel could cause spillover in
    spill_cols <- .get_spill_chs(ms$new, mets$new)
    
    first <- TRUE
    for (i in order(ms$new)) {
        # check if the provided spillovermatrix had the 
        # current channel measured as a spill emitter
        is_new_emitting <- !chs$new[i] %in% chs$emitting
        cur_spillms <- ms$new[spill_cols[[i]]]
        
        if (is_new_emitting) {
            # if channel was never measured, no information is present
            # for spill in any of the potential spill receiver channels
            mass_new_rec <- cur_spillms
        } else {
            # if it has been measured, no information is only present
            # for the channels which were not considered spill receivers 
            # in the spillover matrix provided
            mass_new_rec <- cur_spillms[!cur_spillms %in% ms$receiving]
        }
        if (length(mass_new_rec) > 0) {
            if (first) {
                message("Compensation is likely to be inaccurate.\n",
                    "Spill values for the following interactions\n",
                    "have not been estimated:")
                first <- FALSE
            }
            message(chs$new[i], " -> ", paste(chs$new[
                ms$new %in% mass_new_rec], collapse = ", "))
        }
    }
}