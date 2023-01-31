# ==============================================================================
# get indices of dna columns
# ------------------------------------------------------------------------------
.get_dna_cols <- function(chs, dna) {
    ms <- .get_ms_from_chs(chs)
    n_dna <- length(dna)
    dna_cols <- which(ms %in% dna)
    if (length(dna_cols) != n_dna)
        stop("Not all dna channels found.")
    return(dna_cols)
}

# ==============================================================================
# get indices of bead columns
# ------------------------------------------------------------------------------
.get_bead_cols <- function(chs, beads) {
    ms <- .get_ms_from_chs(chs)
    bead_ms <- if (is.character(beads)) {
        switch(beads, 
            dvs = c(140, 151, 153, 165, 175),
            beta = c(139, 141, 159, 169, 175))
    } else beads
    n_beads <- length(bead_ms)
    bead_cols <- which(ms %in% bead_ms)
    if (length(bead_cols) != n_beads)
        stop("Not all bead channels found.")
    return(bead_cols)
}

# ==============================================================================
# get beads and remove bead-bead doublets
# ------------------------------------------------------------------------------
.get_bead_inds <- function(x, y, assay = "exprs") {
    x <- assignPrelim(x, y, assay = assay, verbose = FALSE)
    applyCutoffs(estCutoffs(x), assay = assay)$bc_id == "is_bead"
}

# ==============================================================================
# bead vs. dna scatter
# ------------------------------------------------------------------------------
#' @import ggplot2
#' @importFrom matrixStats colMins colMaxs
#' @importFrom SummarizedExperiment assay rowData
.plot_bead_scatter <- function(x, dna_chs, bead_chs, assay) {
    # downsample to at most 10k events
    x <- x[, sample(ncol(x), min(1e4, ncol(x)))]
    p <- plotScatter(x, 
        chs = c(dna_chs[1], bead_chs), assay = assay, 
        color_by = "is_bead", label = "channel")
    p$facet$params$nrow <- 1
    p$layers[[1]]$aes_params$alpha <- 0.4
    p$layers[[1]]$aes_params$size <- 0.1
    # get bead intensity boundaries
    y <- t(assay(x, assay))[x$is_bead, ]
    colnames(y) <- channels(x)
    gate <- data.frame(
        variable = factor(bead_chs),
        xmin = min(y[, dna_chs[1]]),
        xmax = max(y[, dna_chs[1]]),
        ymin = colMins(y[, bead_chs]),
        ymax = colMaxs(y[, bead_chs]))
    gate[, -1] <- t(t(gate[, -1]) + rep(0.2, 4)*c(-1, 1))
    # round to nearest 0.25 & expand axes by 0.5 in all directions
    maxs <- vapply(as.list(p$data[c(dna_chs[1], "value")]), 
        function(u) ceiling(max(u)/0.25)*0.25+0.5, numeric(1))
    p + scale_color_manual(values = c("darkgrey", "royalblue")) +
        scale_x_continuous(expand = c(0, 0), 
            breaks = seq(0, maxs[1], 2), limits = c(-0.5, maxs[1])) +
        scale_y_continuous(expand = c(0, 0), 
            breaks = seq(0, maxs[2], 2), limits = c(-0.5, maxs[2])) +
        geom_rect(data = gate, inherit.aes = FALSE, 
            col = "blue", fill = NA, linewidth = 0.4,
            aes_string(
                xmin = "xmin", xmax = "xmax", 
                ymin = "ymin", ymax = "ymax")) + 
        coord_flip() + theme(
            legend.position = "none", 
            panel.spacing = unit(0, "mm"))
}

# ==============================================================================
# plot bead intensitites smoothed by conversion to local medians
# ------------------------------------------------------------------------------
#' @import ggplot2
#' @importFrom dplyr bind_rows group_by mutate_at summarise_if
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
.plot_smooth_beads <- function(b, a, t) {
    # downsample when n > 10k
    i <- TRUE
    if ((n <- length(t)) > 1e4)
        i <- sort(sample(n, 1e4))
    l <- lapply(
        list(before = b, after = a), 
        function(u) data.frame(t(u), t)[i, ])
    df <- bind_rows(l, .id = "id")
    # get baselines
    n_beads <- length(beads <- rownames(a))
    bl <- t(vapply(split.data.frame(df, df$id), function(u)
        colMeans(u[, beads]), numeric(n_beads)))
    bl <- data.frame(bl, id = rownames(bl))
    # reformatt & plot
    bl <- melt(bl, id.vars = "id") 
    df <- melt(df, id.vars = c("t", "id"))
    bl$id <- factor(bl$id, levels = names(l))
    df$id <- factor(df$id, levels = names(l))
    ggplot(df, aes_string("t", "value", col = "variable")) +
        facet_wrap("id", ncol = 1) + geom_line(linewidth = 0.8) + 
        geom_hline(data = bl, linewidth = 0.4, lty = 2, show.legend = FALSE,
            aes_string(yintercept = "value", col = "variable")) +
        scale_color_manual(NULL, values = brewer.pal(9, "Set1")) +
        scale_x_continuous(
            expression("time ("*10^6~"ms)"),
            labels = function(u) u/1e6) +
        labs(y = "smoothed intensity") +
        theme_classic() + theme(
            legend.key.height = unit(0, "mm"),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0, "mm"),
            axis.text = element_text(color = "black"),
            strip.background = element_rect(fill = "white"))
}
