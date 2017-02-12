# ==============================================================================
# Compute deltas
# called by 'assignPrelim()'
# ------------------------------------------------------------------------------

getDeltas <- function(data, bc_key, verbose) {
    
    N <- nrow(data)
    # order barcode intensities within ea. event
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
        
    # NON-CONSTANT NUMBER OF 1'S
    # difference b/w the kth and (k–1)th highest 
    # normalized barcode intensities
    } else {
        
        # find largest barcode separation within ea. event 
        # to assign pos. and neg. barcode values
        diffs <- sapply(1:N, function(x) abs(diff(data[x, bc_orders[x, ]])))
        largest_seps <- apply(diffs, 2, which.max)
        deltas <- sapply(1:N, function(x) diffs[largest_seps[x], x])
    }
    return(deltas)
}
