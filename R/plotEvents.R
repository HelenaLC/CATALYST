#' @rdname plotEvents
#' @title Event plot
#' 
#' @description Plots normalized barcode intensities for a given barcode.
#'
#' @param x 
#'   a \code{\link{dbFrame}}.
#' @param which 
#'   \code{"all"}, numeric or character. Specifies which barcode(s) to plot. 
#'   Valid values are IDs that occur as row names in the \code{bc_key} of the 
#'   supplied \code{\link{dbFrame}}, or 0 for unassigned events. 
#' @param n_events 
#'   numeric. Specifies number of events to plot. Defaults to 100.
#' @param out_path 
#'   character string. If specified, outputs will be generated here.
#' @param name_ext 
#'   character string. If specified, will be appended to the file name. 
#' 
#' @return 
#' Plots intensities normalized by population for each barcode specified
#' by \code{which}: Each event corresponds to the intensities plotted on a 
#' vertical line at a given point along the x-axis. Events are scaled to the 
#' 95\% quantile of the population it has been assigned to. Barcodes with 
#' less than 50 event assignments will be skipped; it is strongly recoomended
#' to remove such populations or reconsider their separation cutoffs.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
#' 
#' @references
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @examples
#' data(sample_ff, sample_key)
#' 
#' # view preliminary assignments
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' plotEvents(x = re, which = "D1", n_events = 1000)
#' 
#' # apply deconvolution parameters
#' re <- estCutoffs(re)
#' re <- applyCutoffs(x = re)
#' plotEvents(x = re, which = "D1", n_events = 500)
#'
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette pdf dev.off
# ------------------------------------------------------------------------------

setMethod(f="plotEvents", 
    signature=signature(x="dbFrame"), 
    definition=function(x, which="all", n_events=100, 
        out_path=NULL, name_ext=NULL) {
        
        # check validity of function arguments
        if ("all" %in% which & length(which) > 1) {
            warning("'which' must either be \"all\"", 
                "or a numeric or character\n",
                "corresponding to valid barcode IDs; ",
                "using default value \"all\".")
        } else if (!"all" %in% which) {
            which <- .check_validity_which(
                which, rownames(bc_key(x)), "events")
        }
        if (!is.numeric(n_events) || n_events == 0 || length(n_events) > 1) {
            warning("'n_events' must be a numeric greater than 0",
                " and of length one;\n  using default value 100.")
            n_events <- 100
        }
        
        # get barcode labels: 
        # channel name if barcodes are single-positive, 
        # sample ID and binary code otherwise
        n_chs <- ncol(bc_key(x))
        ids <- sort(unique(bc_ids(x)))
        if ("all" %in% which) 
            which <- ids
        bc_labs <- .get_bc_labs(x)
        if ("0" %in% ids) 
            bc_labs <- c("Unassigned", bc_labs)
        
        # get colors for plotting:
        # interpolate if more than 11 barcodes, 
        # use colors spaced equally along palette elsewise
        pal <- RColorBrewer::brewer.pal(11,"Spectral")
        if (n_chs > 11) {
            cols <- colorRampPalette(pal)(ncol(normed_bcs(x)))
        } else {
            cols <- pal[ceiling(seq(1, 11, length=ncol(normed_bcs(x))))] 
        }
        
        skipped <- NULL
        p <- list()
        for (id in which) {
            inds <- bc_ids(x) == id
            N <- sum(inds)
            # store IDs with no or insufficient events assigned
            if (N < 50) {
                skipped <- c(skipped, id)
                next
            }
            # subsample events if more than 'n_events' assigned 
            if (N > n_events) {
                inds <- sort(sample(which(inds), n_events))
                n <- n_events
            } else {
                n <- N
            }
            title <- bquote(bold(.(bc_labs[match(id, c("0", rownames(
                bc_key(x))))]))*scriptstyle(" ("*.(N)*" events)"))
            p[[length(p) + 1]] <- .plot_events(x, inds, n, cols, title)
        }
        
        # throw warning about populations with less than 50 event assignments
        if (!is.null(skipped))
            warning("Less than 50 events assigned to Barcode ID(s) ", 
                paste(skipped, collapse=", "), ".")
        
        if (length(p) != 0) {
            if (!is.null(out_path)) {
                pdf(width=10, height=5, file=file.path(out_path, 
                    paste0("event_plot", name_ext, ".pdf")))
                for (i in seq_len(length(p))) plot(p[[i]])
                dev.off()
            } else {
                for (i in seq_len(length(p))) plot(p[[i]])
            }
        }
    })