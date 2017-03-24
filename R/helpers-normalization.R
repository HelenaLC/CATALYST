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
    return(bead_cols)
}

# ==============================================================================
# write FCS of normalized data
# ------------------------------------------------------------------------------

outNormed <- function(ff, normed_es, remove_beads, remove, out_path) {
    if (remove_beads) {
        cells <- new("flowFrame",
            exprs=normed_es[!remove, ],
            parameters=flowCore::parameters(ff),
            description=flowCore::description(ff))
        beads <- new("flowFrame",
            exprs=normed_es[remove, ],
            parameters=flowCore::parameters(ff),
            description=flowCore::description(ff))
        if (is.null(out_path)) {
            flowCore::flowSet(cells, beads)
        } else {
            out_nm <-  file.path(out_path, gsub(".fcs", "", 
                flowCore::description(ff)$GUID, TRUE))
            suppressWarnings(flowCore::write.FCS(cells, 
                paste0(out_nm, "_normalized.fcs")))
            suppressWarnings(flowCore::write.FCS(beads, 
                paste0(out_nm, "_beads.fcs")))
        }
    } else {
        normed <- new("flowFrame",
            exprs=normed_es,
            parameters=flowCore::parameters(ff),
            description=flowCore::description(ff))
        if (is.null(out_path)) {
            normed
        } else {
            suppressWarnings(flowCore::write.FCS(normed, 
                paste0(file.path(out_path, gsub(".fcs", "", 
                    flowCore::description(ff)$GUID, TRUE), 
                    "_normalized.fcs"))))
        }
    }
}

# ==============================================================================
# bead scatters with marginal histograms ontop
# ------------------------------------------------------------------------------
plotBeads <- function(es_t, bead_inds, bead_cols, dna_cols, hist, xlab, gate) {
    
    es <- sinh(es_t)*5
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
    thms <- theme(
        aspect.ratio=1,
        panel.grid.major=element_line(color="grey"),
        panel.grid.minor=element_blank(),
        plot.margin=unit(c(0, .5, .5, .5), "cm"),
        axis.ticks.length=unit(.2, "cm"), 
        axis.ticks=element_line(size=.5),
        axis.title=element_text(size=10), 
        axis.text=element_text(size=8))

    # get axis limits
    x_max <- ceiling(max(es_t[, bead_cols])*2)/2
    y_max <- ceiling(max(es_t[, dna_cols ])*2)/2
    
    p <- list()
    if (hist) {
        for (i in 1:n_beads) {
            med <- median(es[bead_inds, bead_cols[i]])
            min <- min(es[bead_inds, bead_cols[i]])
            max <- max(es[bead_inds, bead_cols[i]])
            med_t <- median(es_t[bead_inds, bead_cols[i]])
            min_t <- min(es_t[bead_inds, bead_cols[i]])
            max_t <- max(es_t[bead_inds, bead_cols[i]])
            p[[length(p)+1]] <- ggplot(data.frame(es_t[bead_inds, ])) +  
                geom_histogram(aes_string(x=chs[bead_cols[i]], y="..ncount.."), 
                    binwidth=.1, size=.025, fill="navy", col="cornflowerblue") +
                annotate("segment", min_t, 0, xend=min_t, yend=.1, col="darkturquoise", size=.75) +
                annotate("segment", med_t, 0, xend=med_t, yend=.05, col="darkturquoise", size=.75) +
                annotate("segment", max_t, 0, xend=max_t, yend=.1, col="darkturquoise", size=.75) +
                annotate("text", .5, .75, hjust=0, size=3, col="darkslateblue", fontface="italic",
                    label=paste0(
                        "Median intensity: ", round(med, 2),
                        "\nLower boundary: ", round(min, 2),
                        "\nUpper boundary: ", round(max, 2))) +
                coord_cartesian(xlim=c(0, x_max), expand=FALSE) +
                scale_x_continuous(labels=function(x) format(x, nsmall=1)) +
                scale_y_continuous(breaks=c(0, .5, 1)) +
                ylab("Normalized counts") + theme_classic()  + thms + 
                theme(axis.title.x=element_text(color="white"), aspect.ratio=.5)
            if (i != 1) p[[length(p)]] <- p[[length(p)]] +
                theme(axis.title.y=element_text(color="white"))
        }
    }
    for (i in 1:n_beads) {
        p[[length(p)+1]] <- ggplot(df, aes_string(col="as.factor(id)",
            x=chs[bead_cols[i]], y=chs[dna_cols[1]])) + theme_bw() + thms +
            scale_colour_manual(values=c("0"="black", "1"="navy")) +
            geom_point(alpha=.25, size=1) + guides(color=FALSE) +
            coord_cartesian(xlim=c(0, x_max), ylim=c(0, y_max), expand=FALSE) +
            scale_x_continuous(labels=function(x) format(x, nsmall=1)) +
            scale_y_continuous(labels=function(x) format(x, nsmall=1)) +
            theme(panel.border=element_rect(colour="black", fill=NA, size=1))
        if (gate) 
            p[[length(p)]] <- p[[length(p)]] + geom_rect(data=gates[i, ], 
                aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), 
                col="darkturquoise", fill=NA, lty=3, size=1, inherit.aes=FALSE) 
        if (!xlab) p[[length(p)]] <- p[[length(p)]] + 
                theme(axis.title.x=element_text(color="white"))
        if (i != 1) p[[length(p)]] <- p[[length(p)]] +
                theme(axis.title.y=element_text(color="white"))
        # if (hist) {
        #     gp1 <- ggplot_gtable(ggplot_build(p[[length(p)]]))
        #     gp2 <- ggplot_gtable(ggplot_build(p[[length(p)-n_beads]]))
        #     maxWidth <- grid::unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
        #     gp1$widths[2:3] <- gp2$widths[2:3] <- maxWidth
        #     p[[length(p)]] <- gp1
        #     p[[length(p)-n_beads]] <- gp2
        # }
    }
    return(p)
}

# ==============================================================================
# plot bead intensitites smoothed by conversion to local medians
# ------------------------------------------------------------------------------
plotSmoothed <- function(df, main) {
    n <- length(chs <- colnames(df))
    baseline <- apply(df[, -1], 2, mean)
    p <- ggplot(df) + ggtitle(main) + theme_classic() + 
        labs(x="Time [s]", y="Local median bead intensity") +
        theme(axis.text=element_text(size=10),
            plot.title=element_text(size=16, face="bold", hjust=.5), 
            axis.title=element_text(size=12, face="bold"),
            panel.grid.major=element_line(color="grey"),
            panel.grid.minor=element_blank())
  
    for (i in 2:n) 
        p <- p + geom_point(size=.25, alpha=.5, 
            aes_string(x=chs[1], y=chs[i], color=as.factor(chs[i]))) + 
        geom_hline(show.legend=FALSE, lty=3, size=.75,
            aes_string(yintercept=baseline[i-1], col=as.factor(chs[i])))
    p + scale_color_manual(name=NULL, 
        values=RColorBrewer::brewer.pal(n-1, "Set1")) +
        guides(colour=guide_legend(override.aes=list(size=3)))
}
