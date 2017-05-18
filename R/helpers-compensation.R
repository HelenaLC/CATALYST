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
make_symetric <- function(x) {
    sm <- diag(ncol(x))
    rownames(sm) <- colnames(sm) <- colnames(x)
    sm[rownames(x), colnames(x)] <- x 
    sm
}

# ==============================================================================
# plot for estTrim()
# ------------------------------------------------------------------------------

plot_estTrim <- function(df, trms, xMin, xMax, yMin, yMax, rect, text) {
    opt <- trms[which.min(text$e)]
    ggplot(df, aes_string(x="t", y="m")) +
        geom_vline(aes(xintercept=opt), lty=2, size=.5) +
        geom_jitter(aes_string(fill="Spiller", group="Receiver"),
            col="navy", height=0, width=diff(trms)[1]/5, size=2, alpha=.3) + 
        geom_rect(fill="aliceblue", inherit.aes=FALSE, data=rect, 
            aes_string(xmin="x1", xmax="x2", ymin="y1", ymax="y2")) + 
        geom_text(size=3, col="blue", vjust=.5, data=text, 
            aes_string(label="e", x="x", y="y")) +
        geom_hline(aes(yintercept=0), lty=2, col="red", size=.5) +
        scale_x_continuous(limits=c(xMin, xMax), 
            expand=c(0,0), breaks=trms, labels=format(trms, 2)) +
        scale_y_continuous(limits=c(yMin, yMax+.6), 
            expand=c(0,0), breaks=c(0, yMin:yMax)) +
        labs(x="Trim value used for estimation of spill values", 
            y="Median counts upon compensation") + 
        theme_classic() + theme(legend.position="none", 
            axis.text=element_text(size=8),
            axis.title=element_text(size=10), 
            panel.grid.major.y=element_blank(),
            panel.grid.major.x=element_line(size=.25, color="grey"),
            panel.grid.minor=element_blank())
}
