# ==============================================================================
# get bead columns
# ------------------------------------------------------------------------------
get_bead_cols <- function(channels, beads) {
    ms <- gsub("[[:alpha:][:punct:]]", "", channels)
    if (is.character(beads)) {
        if (beads == "dvs") {
            bead_ms <- c(140, 151, 153, 165, 175)
        } else if (beads == "beta") {
            bead_ms <- c(139, 141, 159, 169, 175)
        }
    } else {
        bead_ms <- beads
    }
    n_beads <- length(bead_ms)
    if (length(bead_cols <- which(ms %in% bead_ms)) != n_beads)
        stop("Not all bead channels found.")
    return(list(bead_cols, bead_ms))
}

# ==============================================================================
# themes for plotting
# ------------------------------------------------------------------------------
thms <- theme(axis.text=element_text(size=10),
    plot.title=element_text(size=16, face="bold", hjust=.5), 
    axis.title=element_text(size=12, face="bold"),
    panel.grid.major=element_line(color="grey"),
    panel.grid.minor=element_blank())

# ==============================================================================
# bead scatters with marginal histograms ontop
# ------------------------------------------------------------------------------
plotBeads <- function(es_t, bead_inds, bead_cols, dna_cols, hist, xlab, gate) {
    
    chs <- colnames(es_t)
    n_beads <- length(bead_cols)
    
    df <- data.frame(es_t[, c(bead_cols, dna_cols)], id=as.numeric(bead_inds))
    if (gate) {
        gates <- data.frame(t(sapply(bead_cols, function(k) c(
            min(df[bead_inds, chs[k]]), 
            max(df[bead_inds, chs[k]]), 
            min(df[bead_inds, chs[dna_cols[1]]]),
            max(df[bead_inds, chs[dna_cols[1]]])))))
        colnames(gates) <- c("xmin", "xmax", "ymin", "ymax")
    }

    # sample 25'000 events for plotting
    N <- 25e3
    if (nrow(df) > N)
        df <- df[sample(nrow(df), N), ]
    
    # themes for plotting
    thms <- ggplot2::theme(
        aspect.ratio=1,
        panel.grid.major=element_line(color="grey"),
        panel.grid.minor=element_blank(),
        plot.margin=unit(c(0,.5,.5,.5), "cm"),
        axis.ticks.length=unit(.2, "cm"), 
        axis.ticks=element_line(size=.5),
        axis.title=element_text(size=10), 
        axis.text=element_text(size=8))

    # get axis limits
    x_max <- ceiling(max(es_t[, bead_cols])*2)/2
    y_max <- ceiling(max(es_t[, dna_cols ])*2)/2
    
    p <- list()
    for (i in 1:n_beads) {
        # bead versus DNA scatter
        p[[length(p)+1]] <- ggplot(df, aes_string(col="as.factor(id)",
            x=chs[bead_cols[i]], y=chs[dna_cols[1]])) +
            scale_colour_manual(values=c("0"="black", "1"="navy")) +
            geom_point(alpha=.25, size=1) + guides(color=FALSE) +
            coord_cartesian(xlim=c(0, x_max), ylim=c(0, y_max), expand=FALSE) +
            scale_x_continuous(labels=function(x) format(x, nsmall=1)) +
            scale_y_continuous(labels=function(x) format(x, nsmall=1)) +
            theme(panel.border=element_rect(colour="black", fill=NA, size=1)) +
            theme_bw() + thms
        if (gate) 
            p[[length(p)]] <- p[[length(p)]] + geom_rect(data=gates[i, ], 
                aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), 
                col="blue", fill=NA, lty=3, size=1, inherit.aes=FALSE) 

        if (!xlab)  p[[length(p)]] <- p[[length(p)]] + labs(x=NULL)
        if (i != 1) p[[length(p)]] <- p[[length(p)]] + labs(y=NULL)
        # histogram of bead counts
        if (hist) {
            p[[length(p)+1]] <- ggplot(data.frame(es_t[bead_inds, ])) +  
                geom_histogram(aes_string(x=chs[bead_cols[i]], y="..ncount.."), 
                    binwidth=.05, size=.02, fill="navy", col="cornflowerblue") +
                coord_cartesian(xlim=c(0, x_max), expand=FALSE) +
                scale_x_continuous(labels=function(x) format(x, nsmall=1)) +
                scale_y_continuous(breaks=c(0, .5, 1)) +
                labs(x=NULL, y="Normalized counts") + 
                theme_classic()  + thms + theme(aspect.ratio=.5)
            if (i != 1) p[[length(p)]] <- p[[length(p)]] + labs(y=NULL)
            gp1 <- ggplot_gtable(ggplot_build(p[[length(p)]]))
            gp2 <- ggplot_gtable(ggplot_build(p[[length(p)-1]]))
            maxWidth <- grid::unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
            gp1$widths[2:3] <- gp2$widths[2:3] <- maxWidth
            p[[length(p)]] <- gp1
            p[[length(p)-1]] <- gp2
        }
    }
    return(p)
}

# ==============================================================================
# plot raw bead intensitites over time
# ------------------------------------------------------------------------------
# plotBeadsVsTime <- function(es, bead_inds, time_col, bead_cols) {
#     df <- data.frame(es[bead_inds, c(time_col, bead_cols)])
#     n <- length(chs <- colnames(df))
#     p <- ggplot(df) + theme_classic() + thms +
#         ggtitle("Beads prior to smoothing") +
#         labs(x="Time [s]", y="Bead intensity") 
#     for (i in 2:n) 
#         p <- p + geom_point(size=.75, alpha=.5, 
#             aes_string(x=chs[1], y=chs[i], color=as.factor(chs[i]))) 
#     p + scale_color_manual(name=NULL, 
#         values=RColorBrewer::brewer.pal(n-1, "Set1")) +
#         guides(colour=guide_legend(override.aes=list(size=3)))
# }

# ==============================================================================
# plot bead intensitites smoothed by conversion to local medians
# ------------------------------------------------------------------------------
plotSmoothedBeads <- function(df, main) {
    n <- length(chs <- colnames(df))
    baseline <- apply(df[, -1], 2, mean)
    p <- ggplot(df) + theme_classic() + thms + 
        ggtitle(main) +
        labs(x="Time [s]", y="Local median bead intensity") 
    for (i in 2:n) 
        p <- p + geom_point(size=.25, alpha=.5, 
            aes_string(x=chs[1], y=chs[i], color=as.factor(chs[i]))) + 
        geom_hline(show.legend=FALSE, lty=3, size=.75,
            aes_string(yintercept=baseline[i-1], col=as.factor(chs[i])))
    p + scale_color_manual(name=NULL, 
        values=RColorBrewer::brewer.pal(n-1, "Set1")) +
        guides(colour=guide_legend(override.aes=list(size=3)))
}
