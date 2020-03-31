# ==============================================================================
# return barcode IDs 
# called by 'assignPrelim()'
# ------------------------------------------------------------------------------
.get_ids <- function(bcs, bc_key, ids, verbose) {
    
    cutoff <- 0 # used to prevent large neg. values from appearing
    # to have sufficient separation from values near zero
    
    N <- ncol(bcs)
    # order barcode intensities within ea. event
    if (verbose) message(" o ordering")
    bc_orders <- apply(bcs, 2, order, decreasing = TRUE)

    # DOUBLET-FILTERING
    # look at k highest and (n-k)-lowest barcode channels
    if (length(unique(rowSums(bc_key))) == 1) { 
        # number of expected pos. barcode intensities
        n_pos_bcs <- sum(bc_key[1, ])    
        
        # get lowest pos. and highest neg. barcode for ea. event
        lowest_pos  <- bcs[cbind(bc_orders[n_pos_bcs, ],   seq_len(N))]
        highest_neg <- bcs[cbind(bc_orders[n_pos_bcs+1, ], seq_len(N))]
        
        if (verbose) message(" o classifying events")
        # assign binary barcode to ea. event
        codes <- apply(bcs, 1, function(x) as.numeric(x >= lowest_pos))

        # assign barcode ID to ea. event
        lookup <- rowSums(2 ^ col(bc_key) * bc_key)
        preids <- rowSums(2 ^ col(codes)  * codes)
        bc_ids <- ids[match(preids, lookup)]
        
        # exclude events whose pos. barcodes are still very low 
        ex <- bcs[cbind(bc_orders[n_pos_bcs, ], seq_len(N))] < cutoff
        bc_ids[is.na(bc_ids) | ex] <- 0
        
        # NON-CONSTANT NUMBER OF 1'S
        # difference b/w the kth and (k–1)th highest 
        # normalized barcode intensities
    } else {
        # find largest barcode separation within ea. event 
        # to assign pos. and neg. barcode values
        diffs <- vapply(seq_len(N), function(i) 
            abs(diff(bcs[bc_orders[, i], i])), 
            numeric(ncol(bc_key)-1))
        largest_seps <- apply(diffs, 2, function(u)
            order(u, decreasing = TRUE)[1])
        pos <- lapply(seq_len(N), function(i) 
            bc_orders[seq_len(largest_seps[i]), i])
        
        if (verbose) message(" o classifying events")
        # assign binary barcode to ea. event
        codes <- t(vapply(seq_len(N), function(i) 
            as.numeric(seq_len(ncol(bc_key)) %in% pos[[i]]),
            numeric(ncol(bc_key))))
        
        # assign barcode ID to ea. event
        lookup <- rowSums(2 ^ col(bc_key) * bc_key)
        preids <- rowSums(2 ^ col(codes)  * codes)
        bc_ids <- ids[match(preids, lookup)]
        
        # exclude events whose pos. barcodes are still very low 
        # (using bcs, not normalized bcs)
        pos_bcs <- lapply(seq_len(N), function(i) bcs[pos[[i]], i])
        ex <- lapply(pos_bcs, function(x) any(x < cutoff))
        bc_ids[is.na(bc_ids) | unlist(ex)] <- 0
    }
    bc_ids
}

# ==============================================================================
# compute barcode separation for ea. event
# called by 'assignPrelim()'
# ------------------------------------------------------------------------------
.get_deltas <- function(x, bc_key, verbose) {
    # order barcode intensities within ea. event
    bc_orders <- apply(x, 2, order, na.last = TRUE, decreasing = TRUE)
    
    # DOUBLET-FILTERING
    # look at k highest and (n-k)-lowest barcode channels
    n <- ncol(x)
    if (length(unique(rowSums(bc_key))) == 1) { 
        
        # number of expected pos. barcode intensities
        n_pos_bcs <- sum(bc_key[1, ])    
        
        # get lowest pos. and highest neg. barcode for ea. event
        lowest_pos  <- x[cbind(bc_orders[n_pos_bcs, ],   seq_len(n))]
        highest_neg <- x[cbind(bc_orders[n_pos_bcs+1, ], seq_len(n))]
        
        # compute separation b/w pos. and neg. barcodes for ea. event
        lowest_pos - highest_neg
        
        # NON-CONSTANT NUMBER OF 1'S
        # difference b/w the kth and (k–1)th highest 
        # normalized barcode intensities
    } else {
        # find largest barcode separation within ea. event 
        # to assign pos. and neg. barcode values
        diffs <- vapply(seq_len(n), function(i) 
            abs(diff(x[bc_orders[, i], i])),
            numeric(ncol(bc_key) - 1))
        largest_seps <- apply(diffs, 2, function(u)
            order(u, decreasing = TRUE)[1])
        vapply(seq_len(n), function(i) 
            diffs[largest_seps[i], i],
            numeric(1))
    }
}