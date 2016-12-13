# ==============================================================================
# Compute debarcoding
# called by 'assignPrelim()'
# ------------------------------------------------------------------------------

db <- function(data, bcs, bc_key, ids, verbose) {
    
    cutoff <- 0 # used to prevent large neg. values from appearing
    # to have sufficient separation from values near zero
    
    N <- nrow(bcs)
    # order barcode intensities within ea. event
    if (verbose) message(" o ordering")
    bc_orders <- t(apply(data, 1, order, decreasing=TRUE))
    
    # DOUBLET-FILTERING
    # look at k highest and (n-k)-lowest barcode channels
    if (length(unique(rowSums(bc_key))) == 1) { 
        
        # number of expected pos. barcode intensities
        n_pos_bcs <- sum(bc_key[1, ])    
        
        # get lowest pos. and highest neg. barcode for ea. event
        lowest_pos  <- data[cbind(1:N, bc_orders[, n_pos_bcs])]
        highest_neg <- data[cbind(1:N, bc_orders[, n_pos_bcs+1])]
        
        # compute separation b/w pos. and neg. barcodes for ea. event
        deltas <- lowest_pos - highest_neg
        
        if (verbose) message(" o classifying events")
        # assign binary barcode to ea. event
        codes <- apply(data, 2, function(x) as.numeric(x >= lowest_pos))
        
        # assign barcode ID to ea. event
        lookup <- rowSums(2 ^ col(bc_key) * bc_key)
        preids <- rowSums(2 ^ col(codes)  * codes)
        bc_ids <- ids[match(preids, lookup)]
        bc_ids[is.na(bc_ids)] <- 0
        
        # exclude events whose pos. barcodes are still very low 
        # (using bcs, not normalized bcs)
        bc_ids[is.na(bc_ids) | 
            bcs[cbind(1:N, bc_orders[, n_pos_bcs])] < cutoff] <- 0
        
        # NON-CONSTANT NUMBER OF 1'S
        # difference b/w the kth and (k–1)th highest 
        # normalized barcode intensities
    } else {
        
        # find largest barcode separation within ea. event 
        # to assign pos. and neg. barcode values
        diffs <- sapply(1:N, function(x) abs(diff(data[x, bc_orders[x, ]])))
        largest_seps <- apply(diffs, 2, which.max)
        deltas <- diffs[largest_seps]
        pos <- sapply(1:N, function(x) bc_orders[x, 1:largest_seps[x]])
        
        if (verbose) message(" o classifying events")
        # assign binary barcode to ea. event
        codes <- t(sapply(1:N, function(x) 
            as.numeric(1:ncol(bc_key) %in% pos[[x]])))
        
        # assign barcode ID to ea. event
        lookup <- rowSums(2 ^ col(bc_key) * bc_key)
        preids <- rowSums(2 ^ col(codes)  * codes)
        bc_ids <- ids[match(preids, lookup)]
        
        # exclude events whose pos. barcodes are still very low 
        # (using bcs, not normalized bcs)
        pos_bcs <- sapply(1:N, function(x) bcs[x, pos[[x]]])
        ex <- lapply(pos_bcs, function(x) any(x < cutoff))
        bc_ids[is.na(bc_ids) | unlist(ex)] <- 0
    }
    return(list(bc_ids=bc_ids, deltas=deltas))
}
