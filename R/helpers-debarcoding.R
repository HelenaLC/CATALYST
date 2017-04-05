# ==============================================================================
# return barcode IDs 
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
        lowest_pos  <- bcs[cbind(seq_len(N), bc_orders[, n_pos_bcs])]
        highest_neg <- bcs[cbind(seq_len(N), bc_orders[, n_pos_bcs+1])]
        
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
        ex <- bcs[cbind(seq_len(N), bc_orders[, n_pos_bcs])] < cutoff
        bc_ids[is.na(bc_ids) || ex] <- 0
        
        # NON-CONSTANT NUMBER OF 1'S
        # difference b/w the kth and (k–1)th highest 
        # normalized barcode intensities
    } else {
        
        # find largest barcode separation within ea. event 
        # to assign pos. and neg. barcode values
        diffs <- sapply(seq_len(N), function(x) 
            abs(diff(bcs[x, bc_orders[x, ]])))
        largest_seps <- apply(diffs, 2, which.max)
        pos <- sapply(seq_len(N), function(x) 
            bc_orders[x, 1:largest_seps[x]])
        
        if (verbose) message(" o classifying events")
        # assign binary barcode to ea. event
        codes <- t(sapply(seq_len(N), function(x) 
            as.numeric(seq_len(ncol(bc_key)) %in% pos[[x]])))
        
        # assign barcode ID to ea. event
        lookup <- rowSums(2 ^ col(bc_key) * bc_key)
        preids <- rowSums(2 ^ col(codes)  * codes)
        bc_ids <- ids[match(preids, lookup)]
        
        # exclude events whose pos. barcodes are still very low 
        # (using bcs, not normalized bcs)
        pos_bcs <- sapply(seq_len(N), function(x) bcs[x, pos[[x]]])
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
# plot distribution of barcode separations & 
# yields as a function of separation cutoffs
# ------------------------------------------------------------------------------

plot_yield <- function(id, x, seps, n_bcs, lgd, bc_labs) {
    if (id == 0) {
        pal <- RColorBrewer::brewer.pal(11, "Spectral")
        if (n_bcs > 11) {
            cols <- colorRampPalette(pal)(n_bcs)
        } else {
            cols <- sample(pal, n_bcs)
        }
        df_h <- data.frame(sep=seps, count=colSums(counts(x)))
        max <- max(df_h$count)
        df_l <- data.frame(sep=rep(seps, n_bcs), 
            yield=max*c(t(yields(x))), bc=rep(1:n_bcs, each=101))
        
        p <- ggplot() + theme_classic() + theme(
            panel.grid.major=element_line(color="lightgrey"),
            panel.grid.minor=element_blank()) +
            scale_x_continuous("Barcode separation", seq(0, 1, .1)) +
            scale_y_continuous("Yield upon debarcoding", 
                seq(0, max, length=5), labels=paste0(seq(0, 100, 25), "%"), 
                sec.axis=sec_axis(~.*1, "Event count", 
                    labels=scientific_10)) +
            geom_bar(data=df_h, aes_string(x="sep", y="count"),  
                width=1/101, stat="identity", fill="lightgrey", col="white") +
            geom_line(data=df_l, size=.75, aes_string(
                x="sep", y="yield", col="as.factor(bc)")) +
            guides(colour=guide_legend(override.aes=list(size=1))) +
            scale_colour_manual(values=cols, name=NULL, labels=bc_labs) 
        if (!lgd)
            p <- p + guides(colour=FALSE)
    } else {
        df_h <- data.frame(sep=seps, count=counts(x)[id, ])
        max <- max(df_h$count)
        df_l <- data.frame(sep=seps, yield=yields(x)[id, ])
        df_c <- data.frame(sep=seps, fit=max*predict(
            smooth.spline(df_l$sep, df_l$yield))$y)
        df_l$yield <- df_l$yield*max
        
        p <- ggplot() + theme_classic() + theme(
            panel.grid.major=element_line(color="lightgrey"),
            panel.grid.minor=element_blank()) + 
            scale_x_continuous("Barcode separation", seq(0, 1, .1)) +
            scale_y_continuous("Yield upon debarcoding", 
                seq(0, max, length=5), labels=paste0(seq(0, 100, 25), "%"), 
                sec.axis=sec_axis(~.*1, "Event count", labels=scientific_10)) +
            geom_bar(data=data.frame(sep=seps, count=counts(x)[id, ]), 
                aes_string(x="sep", y="count"), width=1/101, 
                stat="identity", fill="navy", colour="lightblue") +
            geom_point(data=df_l, aes_string(x="sep", y="yield"), size=3, 
                fill="green", col="darkgreen", pch=21, alpha=.5, stroke=1) + 
            geom_line(data=df_c, aes_string(x="sep", y="fit"), col="green4") +
            geom_vline(lty=3, size=.75, col="red", 
                aes_string(xintercept="sep_cutoffs(x)[id]"))
    }
    p
}

# ==============================================================================
# get barcode labels: 
# channel name if barcodes are single-positive, 
# barcode ID and binary code otherwise 
# ------------------------------------------------------------------------------
get_bc_labs <- function(x) {
    if (sum(rowSums(bc_key(x)) == 1) == nrow(bc_key(x))) {
        colnames(normed_bcs(x))
    } else {
        paste0(rownames(bc_key(x)), ": ", 
            apply(bc_key(x), 1, function(i) paste(i, collapse="")))
    }
}

# ==============================================================================
# retrieve legend from ggplot
# ------------------------------------------------------------------------------
get_legend <- function(p) {
    tmp <- ggplot_gtable(ggplot_build(p)) 
    lgd <- which(vapply(tmp$grobs, function(x) x$name, "") == "guide-box") 
    return(tmp$grobs[[lgd]]) 
}

# ==============================================================================
# scientific annotation
# ------------------------------------------------------------------------------
scientific_10 <- function(x)
    parse(text=gsub("e", " %*% 10^", format(x, scientific=TRUE)))