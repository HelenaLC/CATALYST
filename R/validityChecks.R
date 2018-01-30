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
    available_clusterings <- colnames(metadata(x)$cluster_codes)
    if (!as.character(k) %in% available_clusterings) {
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
}
