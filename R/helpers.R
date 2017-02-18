# ==============================================================================
# returns barcode IDs 
# called by 'assignPrelim()'
# ------------------------------------------------------------------------------
get_ids <- function(bcs, bc_key, ids, verbose) {
    
    cutoff <- 0 # used to prevent large neg. values from appearing
                # to have sufficient separation from values near zero
    
    N <- nrow(bcs)
    # order barcode intensities within ea. event
    if (verbose) message(" o ordering")
    bc_orders <- t(apply(bcs, 1, order, decreasing=TRUE))
    
    # DOUBLET-FILTERING
    # look at k highest and (n-k)-lowest barcode channels
    if (length(unique(rowSums(bc_key))) == 1) { 
        
        # number of expected pos. barcode intensities
        n_pos_bcs <- sum(bc_key[1, ])    
        
        # get lowest pos. and highest neg. barcode for ea. event
        lowest_pos  <- bcs[cbind(1:N, bc_orders[, n_pos_bcs])]
        highest_neg <- bcs[cbind(1:N, bc_orders[, n_pos_bcs+1])]
        
        if (verbose) message(" o classifying events")
        # assign binary barcode to ea. event
        codes <- apply(bcs, 2, function(x) as.numeric(x >= lowest_pos))
        
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
        diffs <- sapply(1:N, function(x) abs(diff(bcs[x, bc_orders[x, ]])))
        largest_seps <- apply(diffs, 2, which.max)
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
    bc_ids
}

# ==============================================================================
# compute barcode separation for ea. event
# called by 'assignPrelim()'
# ------------------------------------------------------------------------------
get_deltas <- function(data, bc_key, verbose) {
    
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
    deltas
}

# ==============================================================================
# computes mahalanobis distances given current separation cutoffs
# called by 'assignPrelim()'
# ------------------------------------------------------------------------------
get_mhl_dists <- function(x) {
    # get channel and barcode masses
    nms <- colnames(x@exprs)
    ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
    
    # find which columns correspond to barcode masses
    bc_cols <- which(ms %in% colnames(x@bc_key))
    n_bcs <- length(bc_cols)
    
    ids <- unique(x@bc_ids)
    ids <- sort(ids[which(ids != 0)])
    
    # extract barcode columns from FCS file
    bcs <- x@exprs[, bc_cols]
    
    # compute mahalanobis distances of all events
    # given current separation cutoff
    mhl_dists <- numeric(nrow(x@exprs))
    for (i in ids) {
        inds <- which(x@bc_ids == i)
        ex <- inds[x@deltas[inds] < x@sep_cutoffs[rownames(x@bc_key) == i]]
        inds <- inds[!(inds %in% ex)]
        x@bc_ids[ex] <- 0
        sub  <- bcs[inds, ]
        if (length(sub) != n_bcs)
            if (nrow(sub) > n_bcs)
                mhl_dists[inds] <- stats::mahalanobis(
                    x=sub, center=colMeans(sub), cov=stats::cov(sub))
    }
    x@mhl_dists <- mhl_dists
    x
}

# ==============================================================================
# retrieve legend from ggplot
# ------------------------------------------------------------------------------
get_legend <- function(p) {
    tmp <- ggplot_gtable(ggplot_build(p)) 
    lgd <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    return(tmp$grobs[[lgd]]) 
}

# ==============================================================================
# scientific annotation
# ------------------------------------------------------------------------------
scientific_10 <- function(x)
    parse(text=gsub("e", " %*% 10^", format(x, scientific=TRUE)))

# ==============================================================================
# get spillover columns
# ------------------------------------------------------------------------------
get_spill_cols <- function(ms, mets) {
    
    iso_tbl <- list(
        La=138:139, Pr=141, Nd=c(142:146, 148, 150), 
        Sm=c(144, 147:150, 152, 154), Eu=c(151, 153), 
        Gd=c(152, 154:158, 160), Dy=c(156, 158, 160:164), 
        Tb=159, Er=c(162, 164, 166:168, 170), Ho=165, 
        Yb=c(168, 170:174, 176), Tm=169, Lu=175:176)
    
    spill_cols <- list()
    for (i in seq_along(ms)) {
        if (is.na(ms[i])) next
        p1 <- m1 <- ox <- iso <- NULL
        if ((ms[i] + 1)  %in% ms) p1 <- which(ms == (ms[i] + 1))
        if ((ms[i] - 1)  %in% ms) m1 <- which(ms == (ms[i] - 1)) 
        if ((ms[i] + 16) %in% ms) ox <- which(ms == (ms[i] + 16))
        iso <- iso_tbl[[mets[i]]]
        iso <- which(ms %in% iso[iso != ms[i]])
        spill_cols[[i]] <- unique(c(m1, p1, iso, ox))
    }
    spill_cols
}

# ==============================================================================
# make spillover matrix symmetrical
# ------------------------------------------------------------------------------
make_symetric <- function(SM) {
    sm <- diag(ncol(SM))
    rownames(sm) <- colnames(sm) <- colnames(SM)
    sm[rownames(SM), colnames(SM)] <- SM 
    sm
}