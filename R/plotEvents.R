#' @rdname plotEvents
#' @title Event plot
#' @description Plots normalized barcode intensities for a given barcode.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param which 
#'   \code{"all"}, numeric or character specifying which barcode(s) to plot. 
#'   Valid values are IDs that occur as rownames in the \code{bc_key} slot
#'   of the input SCE's \code{metadata}, or 0 for unassigned events. 
#' @param assay character string specifying which 
#'   assay data slot to use. One of \code{assayNames(x)}.
#' @param n single numeric specifying the number of events to plot.
#' @param out_path character string. If specified, 
#'   events plots for all barcodes specified via \code{which} 
#'   will be written to a single PDF file in this location.
#' @param out_name character strings specifying 
#'   the output's file name when \code{!is.null(out_path)}; 
#'   should be provided without(!) file type extension.
#' 
#' @return a list of \code{ggplot} objects.
#' 
#' @details 
#' Plots intensities normalized by population for each barcode specified
#' by \code{which}: Each event corresponds to the intensities plotted on a 
#' vertical line at a given point along the x-axis. Events are scaled to the 
#' 95\% quantile of the population it has been assigned to. Barcodes with 
#' less than 50 event assignments will be skipped; it is strongly recommended
#' to remove such populations or reconsider their separation cutoffs.
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
#' plotEvents(sce, which = "D1")
#'
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay assayNames colData
#' @export

plotEvents <- function(x, which = "all", assay = "scaled", 
    n = 1e3, out_path = NULL, out_name = "event_plot") {
    
    stopifnot(is(x, "SingleCellExperiment"),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        is.numeric(n), length(n) == 1, n > 0, n %% 1 == 0)
    if (!is.null(out_path)) 
        stopifnot(
            is.character(out_path), length(out_path) == 1, dir.exists(out_path),
            is.character(out_name), length(out_name) == 1)
    
    # retreive IDs to include & barcode channels
    rownames(x) <- chs <- channels(x)
    n_bcs <- ncol(bc_key <- metadata(x)$bc_key)
    .check_which(which, rownames(bc_key), "events")
    names(ids) <- ids <- unique(x$bc_id)
    bc_chs <- chs[rowData(x)$is_bc]

    labs <- apply(bc_key[ids, ], 1, paste, collapse = "")
    labs <- paste(ids, labs, sep = ": ")
    names(labs) <- ids
    labs["0"] <- "unassigned"
    
    if (isTRUE(which == "all")) which <- ids
    m <- match(c("0", rownames(bc_key)), which, nomatch = 0)
    names(which) <- which <- which[m]
    
    pal <- brewer.pal(11, "Spectral")
    if (n_bcs > 11) {
        pal <- colorRampPalette(pal)(n_bcs)
    } else {
        idx <- seq(1, 11, length.out = n_bcs)
        pal <- pal[ceiling(idx)]
    }
    
    bc_es <- assay(x, assay)[bc_chs, ]
    cs <- split(seq_len(ncol(x)), x$bc_id)
    ns <- vapply(cs, length, numeric(1))
    cs <- lapply(cs, function(u) sample(u, min(length(u), n)))
    ps <- lapply(which, function(id) {
        if (length(cs[[id]]) == 0) return(NULL)
        df <- data.frame(t(bc_es[, cs[[id]]]), check.names = FALSE)
        df$i <- seq_len(nrow(df))
        gg_df <- melt(df, id.vars = "i")
        ggplot(gg_df, aes_string(x = "i", y = "value", col = "variable")) +
            geom_point() + scale_color_manual(NULL, values = pal) +
            scale_x_continuous(limits = c(0, nrow(df) + 1), expand = c(0, 0)) +
            ggtitle(bquote(bold(.(labs[id]))*" ("*.(ns[id])*" events)")) +
            labs(x = NULL, y = NULL) + theme_bw() + theme(
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                axis.text = element_text(color = "black"),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.key.height = unit(2, "mm"))
    })
    ps <- ps[!vapply(ps, is.null, logical(1))]
    if (length(ps) == 0)
        stop("No or insufficient number of events",
            " assigned to the specified barcode population(s).")
    if (!is.null(out_path)) {
        fn <- paste0(out_name, ".pdf")
        fn <- file.path(out_path, fn)
        pdf(fn, width = 8, height = 4)
        lapply(ps, plot); dev.off()
    } else {
        if (length(ps) == 1) ps[[1]] else ps
    }
}