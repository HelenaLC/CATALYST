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
#' # construct SCE & apply arcsinh-transformation
#' data(sample_ff, sample_key)
#' sce <- fcs2sce(x = sample_ff)
#' 
#' # deconvolute samples & estimate separation cutoffs
#' sce <- assignPrelim(x = sce, bc_key = sample_key)
#' sce <- estCutoffs(x = sce)
#' 
#' # plot for specific sample
#' plotYields(x = sce, which = "C1")
#' 
#' # all barcodes summary plot
#' plotYields(x = sce, which = 0, plotly = TRUE)
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#'
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#' @importFrom htmltools save_html
#' @importFrom methods is
#' @importFrom matrixStats rowMaxs
#' @importFrom plotly config ggplotly hide_legend
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom scales scientific
#' @importFrom S4Vectors metadata
#' @export

plotYields <- function(x, which = 0, 
    out_path = NULL, name_ext = NULL, plotly = FALSE) {
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"),
        !is.null(x$bc_id), !is.null(x$delta),
        !is.null(metadata(x)$bc_key),
        is.null(out_path) || (is.character(out_path) 
            & length(out_path) == 1 & dir.exists(out_path)),
        is.null(name_ext) || (is.character(name_ext) & length(name_ext) == 1),
        is.logical(plotly), length(plotly) == 1)
    n_bcs <- length(ids <- rownames(bc_key <- metadata(x)$bc_key))
    which <- .check_validity_which(which, ids, "yields")
    m <- match(c("0", rownames(bc_key)), which, nomatch = 0)
    names(which) <- which <- which[m]
    
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
        numeric(n_seps))
    counts <- vapply(ids, function(id)
        colSums(counts[cs[[id]], ]),
        numeric(n_seps))

    thm <- function(max, ...) {
        list(theme_classic(), theme(...,
            panel.grid.major = element_line(),
            axis.text = element_text(color = "black"),
            legend.key.height = unit(2, "mm")),
            guides(color = guide_legend(override.aes = list(size = 2))),
            scale_x_continuous(
                name = "Barcode separation", 
                limits = c(-0.025, 1.025),
                breaks = seq(0, 1, 0.1),
                expand = c(0, 0)),
            scale_y_continuous(
                name = "Yield upon debarcoding", 
                limits = c(0, max), 
                breaks = seq(0, max, length = 5), 
                labels = paste0(seq(0, 100, 25), "%"),
                sec.axis = sec_axis(~.*1, "Cell count", labels = scientific)))
    }
    
    cuts <- metadata(x)$sep_cutoffs
    labs <- apply(bc_key, 1, paste, collapse = "")
    labs <- paste(ids, labs, sep = ": ")
    names(labs) <- ids
    
    ps <- lapply(which, function(id) {
        if (id == "0") {
            h <- data.frame(count = rowSums(counts), cutoff = seps)
            l <- melt(yields * (max <- max(h$count)))
            names(l) <- c("cutoff", "bc_id", "yield")
            if (is.numeric(l$bc_id))
                l$bc_id <- factor(l$bc_id, levels = sort(unique(l$bc_id)))
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
                thm(max, legend.position = "none")
        } else {
            df <- data.frame(cutoff = seps, count = counts[, id]) 
            df$yield <- yields[, id] * (max <- max(df$count))
            p <- ggplot(df, aes_string("cutoff")) +
                geom_bar(aes_string(y = "count"), stat = "identity", 
                    size = 0.2, col = "white", fill = "darkgrey") + 
                geom_line(aes_string(y = "yield"), col = "red") + 
                ggtitle(labs[id]) + thm(max) 
            if (!is.null(cuts[id])) {
                cut <- round(cuts[id], 2)
                y <- format(yields[as.character(cut), id] * 100, digits = 4)
                p <- p + geom_label(
                    label = paste0("cutoff: ", cut, "; yield: ", y, "%"),
                    x = cuts[id] + 0.01, y = max, hjust = 0, col = "blue") +
                    geom_vline(xintercept = cuts[id], col = "blue")
            }
            p
        }
    })
    
    if (plotly) {
        # remove geom_label as it has not
        # been implemented in plotly yet
        ps <- lapply(ps, function(p) {
            geoms <- vapply(p$layers, function(u) 
                class(u$geom)[1], character(1))
            is_geom_label <- grep("label", geoms, ignore.case = TRUE)
            p$layers[is_geom_label] <- NULL
            return(p)
        })
        ps <- lapply(which, function(id) {
            p <- switch(as.character(id), 
                "0" = hide_legend(ggplotly(ps[[id]], tooltip = "text")),
                ggplotly(ps[[id]], tooltip = c("cutoff", "yield", "count")))
            config(p, displayModeBar = FALSE)
        })
    }
    
    if (!is.null(out_path)) {
        ext <- ifelse(plotly, ".html", ".pdf")
        fn <- paste0("yield_plot", name_ext, ext)
        fn <- file.path(out_path, fn)
        if (plotly) {
            save_html(ps, fn)
        } else {
            pdf(fn, width = 7, height = 3.5)
            lapply(ps, print); dev.off()
        }
    } else {
        if (length(ps) == 1) ps[[1]] else ps
    }
}