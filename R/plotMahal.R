#' @rdname plotMahal
#' @title Biaxial plot
#' 
#' @description 
#' Histogram of counts and plot of yields as a function of separation cutoffs.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param which character string. Specifies which barcode to plot.
#' @param assay character string specifying which assay to use.
#' @param n numeric. Number of cells to subsample; use NULL to include all. 
#' 
#' @return Plots all inter-barcode interactions for the population specified
#' by argument \code{which}. Events are colored by their Mahalanobis distance. 
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @examples
#' data(sample_ff, sample_key)
#' sce <- prepData(sample_ff, by_time = FALSE)
#' sce <- assignPrelim(sce, sample_key)
#' sce <- estCutoffs(sce)
#' sce <- applyCutoffs(sce)
#' plotMahal(sce, which = "B3")
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange 
#' @importFrom cowplot get_legend
#' @importFrom methods is
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SummarizedExperiment assayNames
#' @export

plotMahal <- function(x, which, assay = "exprs", n = 1e3) {
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        is.null(n) || (is.numeric(n) & length(n) == 1 & n > 1))
    n_bcs <- ncol(bc_key <- metadata(x)$bc_key)
    if (ncol(bc_key) > 7)
        stop("\nToo many barcodes: ",
            "Plotting all inter-barcode interactions is infeasible.\n",
            "Using the default 'mhl_cutoff' value of 30 is recommended.")
    stopifnot(is.character(which), length(which) == 1, which %in% rownames(bc_key))

    # subset specified barcode population 
    # & (optionally) subsample cells
    cs <- which(x$bc_id == which & !is.na(x$mhl_dist))
    if (!is.null(n) && n < length(cs))
        cs <- sample(cs, min(n, length(cs)))
    ds <- x$mhl_dist[cs]
    chs <- rownames(x)
    ms <- .get_ms_from_chs(chs)
    m <- match(colnames(bc_key), ms)
    bc_chs <- rownames(x)[m]
    y <- assay(x, assay)[bc_chs, cs]
    df <- data.frame(t(y), d = ds, check.names = FALSE)
    
    # specify plotting aesthetics
    thm <- theme_classic() + theme(
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1)
    
    # get axes limits
    max <- ceiling(max(ds, na.rm = TRUE) / 5) * 5
    axes_min <- floor(min(y, na.rm = TRUE) / 0.1) * 0.1
    axes_max <- ceiling(max(y, na.rm = TRUE) / 0.1) * 0.1
    lims <- c(axes_min, axes_max)
    
    # initialize plot list
    n_ps <- sum(seq_len(n_bcs))
    ps <- vector("list", n_ps)
    
    # histogram of barcode distributions
    hi <- c(1, cumsum(seq(n_bcs, 2)) + 1)
    for (p in hi) {
        ch <- sprintf("`%s`", bc_chs[match(p, hi)])
        ps[[p]] <- ggplot(df) + 
            scale_x_continuous(limits = lims) +
            geom_histogram(aes_string(x = ch), 
                breaks = seq(axes_min + 0.05, axes_max - 0.05, 0.1),
                binwidth = 0.1, fill = "black", color = NA) + 
            labs(x = " ", y = " ") + coord_fixed(1) + thm
    }
    m <- matrix(NA, n_bcs, n_bcs)
    m[lower.tri(m, diag = TRUE)] <- seq_along(ps)
    
    # bead vs. bead scatters color coded 
    # according to Mahalanobis distances
    first <- TRUE
    for (i in seq_len(n_bcs - 1)) {
        for (j in seq((i + 1), n_bcs)) {
            chs <- sprintf("`%s`", bc_chs[c(i, j)])
            ps[[m[j, i]]] <- ggplot(df) + geom_point(size = 0.8,
                    aes_string(x = chs[1], y = chs[2], col = "d")) +
                guides(color = "none") + scale_color_gradientn(
                    sprintf("%s: %s", which, 
                        paste(bc_key[which, ], collapse = "")),
                    colors = rev(brewer.pal(11, "RdYlBu")),
                    limits = c(0, max), breaks = seq(0, max, 5)) +
                lims(x = lims, y = lims) + labs(x = " ", y = " ") + thm
            if (first) {
                # get color bar
                foo <- ps[[m[j, i]]] + guides(colour = guide_colourbar(
                    title.position = "top", title.hjust = 0.5)) + 
                    theme(legend.direction = "horizontal",
                        legend.title = element_text(face="bold"),
                        legend.key.height = unit(0.5, "line"),
                        legend.key.width = unit(4, "line"),
                        legend.text = element_text(size = 8), 
                        legend.key = element_blank())
                suppressWarnings(lgd <- get_legend(foo))
                first <- FALSE
            }
        }
    }
    # add axes labels for left column & bottom row plots
    for (i in seq_len(n_bcs)) {
        l <- m[i, 1]; b <- m[n_bcs, i]
        ps[[b]] <- ps[[b]] + labs(x = bc_chs[i])
        ps[[l]] <- ps[[l]] + labs(y = bc_chs[i])
    }
    ps[[n_ps + 1]] <- lgd

    m <- rbind(rep(n_ps + 1, n_bcs), m)
    hs <- c(2.5, rep(5, n_bcs)); ws  <- rep(5, n_bcs)
    grid.arrange(grobs = ps, layout_matrix = m, heights = hs, widths = ws)
}
