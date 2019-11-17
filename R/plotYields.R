#' @rdname plotYields
#' @title Yield plot
#' @description Distribution of barcode separations 
#' and yields as a function of separation cutoffs.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param which 0, numeric or character. Specifies which barcode(s) to plot. 
#'   Valid values are IDs that occur as row names of \code{bc_key(x)}; 
#'   0 (the default) will generate a summary plot with all barcodes.
#' @param out_path character string. If specified, outputs will be generated here.
#' @param name_ext character string. If specified, will be appended to the plot's name. 
#' @param plotly logical. Should an interactive plot be rendered?
#' 
#' @return plots the distribution of barcode separations and yields upon 
#' debarcoding as a function of separation cutoffs. If available, currently 
#' used separation cutoffs as well as their resulting yields will be indicated 
#' in the plot`s main title.
#' 
#' @details
#' The overall yield that will be achieved upon application of the specified 
#' set of separation cutoffs is indicated in the summary plot. Respective 
#' separation thresholds and their resulting yields are included in each 
#' barcode's plot. The separation cutoff value should be chosen such that
#' it appropriately balances confidence in barcode assignment and cell yield.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' re <- estCutoffs(x = re)
#' 
#' # all barcodes summary plot
#' plotYields(x = re, which = 0)
#' 
#' # plot for specific sample
#' plotYields(x = re, which = "C1")
#' 
#' @author Helena L. Crowell
#'
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @importFrom ggplot2
#' @importFrom htmltools save_html
#' @importFrom methods is
#' @importFrom matrixStats rowMaxs
#' @importFrom plotly config ggplotly hide_legend
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment metadata
#' @export

plotYields <- function(x, which = 0, 
    out_path = NULL, name_ext = NULL, plotly = FALSE) {
    stopifnot(is(x, "SingleCellExperiment"))
    
    bc_key <- metadata(x)$bc_key
    ids <- rownames(bc_key)
    n_bcs <- nrow(bc_key)
    
    # compute yields & cell counts
    n_seps <- length(names(seps) <- seps <- seq(0, 1, 0.01))
    yields <- vapply(seps, function(u) x$delta >= u, numeric(ncol(x)))
    counts <- vapply(seps, function(u) {
        v <- ifelse(u == seps[n_seps], u, u + 0.01)
        yields[, as.character(u)] & x$delta < v
    }, numeric(ncol(x)))
    
    # split cell by barcode ID
    cs <- split(seq_len(ncol(x)), x$bc_id)
    yields <- vapply(ids, function(id)
        colMeans(yields[cs[[id]], ]),
        numeric(length(seps)))
    counts <- vapply(ids, function(id)
        colSums(counts[cs[[id]], ]),
        numeric(length(seps)))

    thm <- function(max) {
        list(theme_classic(), theme(
            panel.grid.major = element_line(),
            axis.text = element_text(color = "black"),
            legend.key.height = unit(2, "mm")),
            guides(color = guide_legend(override.aes = list(size = 2))),
            scale_x_continuous(
                name = "Barcode separation", 
                limits = c(-0.025, 1.025),
                breaks = seq(0, 1, 0.1)),
            scale_y_continuous(
                name = "Yield upon debarcoding", 
                limits = c(0, max), 
                breaks = seq(0, max, length = 5), 
                labels = paste0(seq(0, 100, 25), "%"),
                sec.axis = sec_axis(~.*1, "Cell count", labels = scientific)))
    }
    
    labs <- apply(bc_key, 1, paste, collapse = "")
    labs <- paste(ids, labs, sep = ": ")
    names(labs) <- ids
    
    ps <- lapply(which, function(id) {
        if (id == "0") {
            h <- data.frame(count = rowSums(counts), cutoff = seps)
            l <- melt(yields * (max <- max(h$count)))
            names(l) <- c("cutoff", "bc_id", "yield")
            pal <- brewer.pal(11, "Spectral")
            if (n_bcs > 11) {
                pal <- colorRampPalette(pal)(n_bcs)
            } else {
                pal <- sample(pal, n_bcs)
            }
            ggplot() +
                geom_bar(data = h, aes_string("cutoff", "count"), 
                    stat = "identity", width = 1 / n_seps, 
                    fill = "lightgrey", col = "white") +
                geom_line(data = l, aes_string("cutoff", "yield", col = "bc_id")) + 
                scale_color_manual(NULL, values = pal, labels = labs) +
                thm(max)
        } else {
            df <- data.frame(cutoff = seps, count = counts[, id]) 
            df$yield <- yields[, id] * (max <- max(df$count))
            ggplot(df, aes_string("cutoff")) +
                geom_bar(aes_string(y = "count"), stat = "identity") + 
                geom_line(aes_string(y = "yield")) + 
                ggtitle(labs[id]) + thm(max) 
        }
    })
    
    if (plotly)
        ps <- lapply(which, function(id)
            p <- switch(as.character(id), 
                "0" = hide_legend(ggplotly(ps[[i]], tooltip = "text")),
                ggplotly(ps[[i]], tooltip = c("cutoff", "yield", "count")))
            config(p, displayModeBar = FALSE))
    
    if (!is.null(out_path)) {
        fn <- paste0("yield_plot", name_ext, ".html")
        save_html(ps, file.path(out_path, fn))
    } else {
        ps
    }
}