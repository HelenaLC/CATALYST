# ==============================================================================
# Validity check for 'which' in 'plotEvents()' and 'plotYields()'
#       - stop if not a single ID is valid
#       - warning if some ID(s) is/are not valid and remove it/them
# ------------------------------------------------------------------------------
check_validity_which <- function(which, ids, fct) {
    
    msg_events <- c(
        " Valid values for 'which' are IDs that occur as row names in the\n",
        " 'bc_key' slot of the supplied 'dbFrame', or 0 for unassigned events.")
    
    msg_yields <- c(
        " Valid values for 'which' are IDs that occur as row names in the\n",
        " 'bc_key' slot of the supplied 'dbFrame', or 0 for all barcodes.")
    
    if (length(which) == 1 && !which %in% c(0, ids)) {
        if (fct == "events") {
            stop(paste(which), 
                " is not a valid barcode ID.\n", 
                msg_events, call.=FALSE)
        } else if (fct == "yields") {
            stop(paste(which), 
                " is not a valid barcode ID.\n", 
                msg_yields, call.=FALSE)
        }
    } else {
        new <- which[!is.na(match(which, c(0, ids)))]
        removed <- which[!which %in% new]
        if (length(new) == 0) {
            if (fct == "events") {
                stop(paste(removed, collapse=", "), 
                    " are not valid barcode IDs.\n",
                    msg_events, call.=FALSE)
            } else if (fct == "yields") {
                stop(paste(removed, collapse=", "), 
                    " are not valid barcode IDs.\n",
                    msg_yields, call.=FALSE)
            }
        } else if (length(new) != length(which)) {
            which <- new
            if (length(removed) == 1) {
                warning(paste(removed),
                    " is not a valid barcode ID and has been skipped.",
                    call.=FALSE)
            } else {
                warning(paste(removed, collapse=", "),
                    " are not valid barcode IDs and have been skipped.",
                    call.=FALSE)
            }
        } 
    }
    as.character(which)
}

# ==============================================================================
# Validity check for clustering 'k': should be...
#       - a numeric that accesses the FlowSOM clustering (100),
#         or ConsensusClusterPlus metaclustering (2, ..., 20)
#       - a character string that matches with a 'label' specifying
#         a merging done with 'mergeClusters'
# ------------------------------------------------------------------------------
check_validity_of_k <- function(x, k) {
    available_clusterings <- colnames(S4Vectors::metadata(x)$cluster_codes)
    is_valid_clustering <- as.character(k) %in% available_clusterings
    if (!is.null(k) && !is_valid_clustering) {
        if (is.numeric(k)) {
            txt <- k 
        } else {
            txt <- dQuote(k)
        }
        ks <- suppressWarnings(as.numeric(colnames(cluster_codes(x))))
        ks <- ks[!is.na(ks)]
        stop("Clustering 'k = ", txt, "' doesnt't exist. ", 
            "Should be one of\n  ", paste(c(ks, dQuote(setdiff(
                available_clusterings, ks))), collapse=", "))
    }
    return(as.character(k))
}

# ==============================================================================
# Validity check for columns (used in daFrame constructor): should be...
#        - a logical vector
#        - numeric vector of indices
#        - character vector of column names
# ------------------------------------------------------------------------------
check_validity_cols <- function(cols_to_use, col_nms) {
    n_cols <- length(col_nms)
    
    check1 <- is.logical(cols_to_use) & length(cols_to_use) == n_cols
    check2 <- all(cols_to_use %in% col_nms)
    check3 <- TRUE
    if (is.integer(cols_to_use))
        check3 <- min(cols_to_use >= 1) && max(cols_to_use <= ncol)

    if (!any(check1, check2, check3))
        stop("Invalid argument 'cols_to_use'. Should be either", 
            " a logial vector,\n  a numeric vector of indices, or",
            " a character vector of column names.")
}
