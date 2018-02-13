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
        bc_ids[is.na(bc_ids) | ex] <- 0
        
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
plot_yields <- function(id, x, seps, n_bcs, bc_labs) {
    if (id == 0) {
        pal <- RColorBrewer::brewer.pal(11, "Spectral")
        if (n_bcs > 11) {
            cols <- colorRampPalette(pal)(n_bcs)
        } else {
            cols <- sample(pal, n_bcs)
        }
        counts <- colSums(counts(x))
        max <- max(counts)
        hist <- data.frame(seps, counts)
        line <- reshape2::melt(data.frame(
            seps, yields=t(yields(x))*max), id.var=seps)
        line$Sample <- rep(bc_labs, each=length(seps))
        
        p <- suppressWarnings(ggplot() + 
            geom_bar(data=hist, aes_string(x="seps", y="counts"), width=1/101, 
                stat="identity", fill="lightgrey", col="white") +
            geom_line(data=line, size=.5, aes_string(x="seps", y="value",
                col="as.factor(variable)", text="Sample")) +
            guides(colour=guide_legend(override.aes=list(size=1))) +
            scale_colour_manual(values=cols, name=NULL, labels=bc_labs))
    } else {
        max <- max(counts(x)[id, ])
        df <- data.frame(Cutoff=seps, 
            Count=counts(x)[id, ], 
            yield=yields(x)[id, ]*max,
            fit=max*predict(smooth.spline(seps, yields(x)[id, ]), seps)$y)
        df$Yield <- paste0(sprintf("%2.2f", round(100*df$yield/max, 2)), "%")
        inds <- seq(1, nrow(df), 2)
        
        p <- ggplot(df) + 
            geom_bar(aes_string(x="Cutoff", y="Count"), width=1/101, size=.3, 
                stat="identity", fill="lavender", colour="darkslateblue") +
            geom_point(data=data.frame(df[inds, ]),
                aes_string(x="Cutoff", y="yield", group="Yield"), 
                fill="mintcream", color="aquamarine4", 
                pch=21, size=3, stroke=.75) + 
            geom_line(data=data.frame(x=seps, y=df$fit), 
                aes_string(x="x", y="y"), 
                col="aquamarine3", size=.75) +
            geom_vline(lty=3, size=.75, col="red2", 
                aes_string(xintercept="sep_cutoffs(x)[id]")) 
    }
    p + scale_x_continuous(name="Barcode separation", 
        breaks=seq(0, 1, .1), limits=c(-.025,1.025), expand=c(0,0)) +
        scale_y_continuous(name="Yield upon debarcoding", 
            breaks=seq(0, max, length=5), labels=paste0(seq(0, 100, 25), "%"),
            limits=c(-.05*max, max+.05*max), expand=c(0,0), sec.axis=sec_axis(
                trans=~.*1, name="Event count", labels=scientific_10)) +
        theme_classic() + theme(axis.text=element_text(color="black"),
            panel.grid.major=element_line(color="grey", size=.25))
}

# ==============================================================================
# helper for plotEvents()
# ------------------------------------------------------------------------------
plot_events <- function(x, inds, n, cols, title) {
    df <- reshape2::melt(data.frame(
        x=seq_len(n), y=normed_bcs(x)[inds, ]), id.var="x")
    y_max <- ceiling(4*max(df$value))/4
    ggplot(df, aes_string(x="x", y="value", col="as.factor(variable)")) +
        geom_point(stroke=.5, size=1+100/n, alpha=.75) +
        scale_colour_manual(values=cols, labels=colnames(normed_bcs(x))) + 
        guides(colour=guide_legend(title=NULL,
            override.aes=list(alpha=1, size=2.5))) +
        scale_x_continuous(limits=c(0, n+1), 
            breaks=seq_len(n), expand=c(0, 0)) + 
        scale_y_continuous(limits=c(-.125, y_max+.125), 
            breaks=seq(0, y_max, .25), expand=c(0, 0)) +
        labs(x=NULL, y="Normalized intensity") + 
        ggtitle(title) + theme_classic() + theme(
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_text(color="black"),
            panel.grid.major.y=element_line(color="grey", size=.25),
            panel.grid.minor=element_blank())
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
    
