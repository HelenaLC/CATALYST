# ==============================================================================
# get indices of bead columns
# ------------------------------------------------------------------------------
.get_bead_cols <- function(channels, beads) {
    ms <- gsub("[[:alpha:][:punct:]]", "", channels)
    if (is.character(beads)) {
        if (isTRUE(beads == "dvs")) {
            bead_ms <- c(140, 151, 153, 165, 175)
        } else if (isTRUE(beads == "beta")) {
            bead_ms <- c(139, 141, 159, 169, 175)
        }
    } else {
        bead_ms <- beads
    }
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
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment assay colData
.plot_bead_scatter <- function(x, dna_chs, bead_chs, assay = "exprs") {
    dna <- dna_chs[1]
    df <- data.frame(
        i = seq_len(ncol(x)), colData(x),
        t(assay(x, assay)[c(bead_chs, dna), ]))
    df <- melt(df, id.vars = c("i", names(colData(x)), dna))
    if(ncol(x) > 1e4) df <- df[df$i %in% sample(ncol(x), 1e4), ]
    ggplot(df, aes_string(x = "value", y = dna, col = "is_bead")) +
        facet_wrap("variable", ncol = length(bead_chs)) + 
        geom_rect(data = df[df$is_bead, ], aes_string(
            xmin = "min(value)", xmax = "max(value)",
            ymin = "min(get(dna))", ymax = "max(get(dna))"),
            fill = NA, size = 0.1, show.legend = FALSE) +
        geom_point(pch = 16, size = 0.5, alpha = 0.5, show.legend = FALSE) + 
        scale_color_manual(NULL, values = c("black", "red")) +
        coord_equal() + theme_bw() + theme(
            panel.grid = element_blank(),
            axis.title.x = element_blank(),
            axis.text = element_text(color = "black"),
            panel.spacing = unit(0, "mm"),
            strip.background = element_rect(fill = "white"))
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
        facet_wrap("id", ncol = 1) + geom_line(size = 0.8) + 
        geom_hline(data = bl, size = 0.4, lty = 2, show.legend = FALSE,
            aes_string(yintercept = "value", col = "variable")) +
        scale_color_manual(NULL, values = brewer.pal(9, "Set1")) +
        labs(x = "time (s)", y = "smoothed intensity") +
        theme_classic() + theme(
            legend.key.height = unit(0, "mm"),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0, "mm"),
            axis.text = element_text(color = "black"),
            strip.background = element_rect(fill = "white"))
}
