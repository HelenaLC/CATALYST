# ==============================================================================
# Plot events
# ------------------------------------------------------------------------------
#' @rdname plotEvents
#' @title Event plot
#' 
#' @description 
#' Shows normalized barcode intensities for a given barcode.
#'
#' @param x a \code{\link{dbFrame}}.
#' @param which 
#' "all", numeric or character. Specifies which barcode(s) to plot. 
#' Valid values are IDs that occur as row names in the \code{bc_key} 
#' of the supplied \code{\link{dbFrame}}, or 0 for unassigned events. 
#' Defaults to "all".
#' @param n_events 
#' numeric. Specifies number of events to plot. Defaults to 100.
#' @param out_path 
#' a character string. If specified, outputs will be generated 
#' in this location. Defaults to NULL.
#' @param name_ext 
#' a character string. If specified, will be appended to the plot's name. 
#' Defaults to NULL.
#' 
#' @return 
#' plots intensities normalized by population for each barcode specified
#' by \code{which}: Each event corresponds to the intensities plotted on a 
#' vertical line at a given point along the x-axis. Events are scaled to the 
#' 95\% quantile of the population it has been assigned to. Barcodes with 
#' less than 50 event assignments will be skipped; it is strongly recoomended
#' to remove such populations or reconsider their separation cutoffs.
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
#' @references
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette pdf dev.off
# ==============================================================================

setMethod(f="plotEvents", 
    signature=signature(x="dbFrame"), 
    definition=function(x, which="all", n_events=100, 
        out_path=NULL, name_ext=NULL) {
        
        # ······································································
        # check validity of function arguments
        if ("all" %in% which & length(which) > 1) {
            warning("'which' must either \"all\", or a numeric or character\n",
                "corresponding to valid barcode IDs; ",
                "using default value \"all\".")
        } else if (!"all" %in% which) {
            which <- check_validity_which(
                which, rownames(bc_key(x)), "events")
        }
        if (!is.numeric(n_events) || n_events == 0 || length(n_events) > 1) {
            warning("'n_events' must be a numeric greater than 0",
                " and of length one;\n  using default value 100.")
            n_events <- 100
        }
        # ······································································
        
        n_chs <- ncol(bc_key(x))
        ids <- sort(unique(bc_ids(x)))
        if ("all" %in% which) 
            which <- ids
        bc_labs <- get_bc_labs(x)
        if (0 %in% ids) 
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
            n <- sum(inds)
            # store IDs with no or insufficient events assigned
            if (n < 50) {
                skipped <- c(skipped, id)
                next
            }
            # subsample events if more than 'n_events' assigned 
            if (n > n_events) 
                inds <- sort(sample(which(inds), n_events))
            # use normalized barcode intensities
            ints <- normed_bcs(x)[inds, ]
            
            df <- data.frame(
                event=rep(seq_len(length(ints) / n_chs), each=n_chs),
                intensity=c(t(ints)),
                bc=rep(seq_len(n_chs), (length(ints) / n_chs)))
            
            p[[length(p) + 1]] <- ggplot(df) + 
                geom_point(stroke=.5, size=1+100/sum(inds), aes_string(
                    x="event", y="intensity", col="as.factor(bc)", alpha=.8)) +
                scale_colour_manual(values=cols, name=NULL, 
                    labels=colnames(normed_bcs(x))) + 
                guides(alpha=FALSE, 
                    colour=guide_legend(override.aes=list(size=3))) +
                scale_x_discrete(limits=NULL, labels=NULL) + 
                expand_limits(x=c(0, nrow(ints)+1)) +
                ylim(floor(4*min(ints))/4, ceiling(4*max(ints))/4) + 
                xlab("Event number") + ylab("Normalized intensity") + 
                ggtitle(bquote(bold(.(bc_labs[match(id, rownames(
                    bc_key(x)))]))*scriptstyle(" ("*.(n)*" events)"))) +
                theme_bw() + theme(legend.key=element_blank(),
                    panel.grid.major=element_line(color="lightgrey"),
                    panel.grid.minor=element_blank())
        }
        
        # ······································································
        # throw warning about populations with less than 50 event assignments
        if (!is.null(skipped))
            warning("Less than 50 events assigned to Barcode ID(s) ", 
                paste(skipped, collapse=", "), ".")
        # ······································································
        
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