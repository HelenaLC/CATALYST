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
#' a character string. If specified, outputs will be generated in this location. 
#' Defaults to NULL.
#' @param name_ext 
#' a character string. If specified, will be appended to the plot's name. 
#' Defaults to NULL.
#' 
#' @details 
#' EXPLAIN NORMALIZATION HERE
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' plotEvents(x = re, which = "B3", n_events = "all")
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' plotEvents(x = re, which = "B3", n_events = "all")
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
#' @export
# ==============================================================================

setMethod(f="plotEvents", 
    signature=signature(x="dbFrame"), 
    definition=function(x, which="all", n_events=100, 
        out_path=NULL, name_ext=NULL) {
        
        # ----------------------------------------------------------------------
        # check validity of function arguments
        if (!"all" %in% which) 
            which <- check_validity_which(
                which, rownames(x@bc_key), "events")
        if (!is.numeric(n_events) || n_events == 0 || length(n_events) > 1) {
            warning("'n_events' must be a numeric greater than 0 and ",
                "of length one;\n using default value 100.", call.=FALSE)
            n_events <- 100
        }
        # ----------------------------------------------------------------------
        
        set.seed(8)
        
        nms <- colnames(x@exprs)
        ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
        
        ids <- sort(unique(x@bc_ids))
        n_bcs <- nrow(x@bc_key)
        n_chs <- ncol(x@bc_key)
        
        if (sum(rowSums(x@bc_key) == 1) == n_bcs) {
            bc_labs <- paste(colnames(x@exprs))[
                ms %in% as.numeric(colnames(x@bc_key))]
        } else {
            bc_labs <- paste0(rownames(x@bc_key), ": ", 
                apply(x@bc_key, 1, function(x) paste(x, collapse="")))
        }
        if (0 %in% ids) 
            bc_labs <- c("Unassigned", bc_labs)
        bc_labs <- bc_labs[c(0, rownames(x@bc_key)) %in% ids]
        
        labs <- colnames(x@normed_bcs)
        pal <- RColorBrewer::brewer.pal(11,"Spectral")
        if (n_chs > 11) {
            cols <- colorRampPalette(pal)(ncol(x@normed_bcs))
        } else {
            cols <- pal[ceiling(seq(1, 11, length=ncol(x@normed_bcs)))] 
        }
        thms <- theme_bw() + theme(legend.key=element_blank(),
            panel.grid.major=element_line(color="lightgrey"),
            panel.grid.minor=element_blank())
        
        if ("all" %in% which) {
            which_ids <- ids
        } else {
            which_ids <- which
        }
        
        # initialize vector to store IDs with no or insufficient events assigned
        skipped <- NULL
        
        p <- list()
        for (id in which_ids) {
            inds <- which(x@bc_ids == id)
            n <- length(inds)
            if (length(inds) < 50) {
                skipped <- c(skipped, id)
                next
            }
            if (length(inds) > n_events) 
                inds <- sort(sample(inds, n_events))
            psize <- 1 + 100 / length(inds)
            sub <- x@normed_bcs[inds, ]
            
            tcks <- c(1,nrow(sub))
            ymin <-   floor(4 * min(sub)) / 4
            ymax <- ceiling(4 * max(sub)) / 4
            
            df <- data.frame(
                event=rep(1:(length(sub) / n_chs), each=n_chs),
                intensity=c(t(sub)),
                bc=rep(1:n_chs, (length(sub) / n_chs)))

            p[[length(p) + 1]] <- ggplot(df) + 
                geom_point(stroke=.5, size=psize, aes_string(
                    x="event", y="intensity", col="as.factor(bc)", alpha=.8)) +
                scale_colour_manual(values=cols, name=NULL, labels=labs) + 
                guides(alpha=FALSE, 
                    colour=guide_legend(override.aes=list(size=3))) +
                scale_x_discrete(limits=NULL, labels=NULL) + ylim(ymin, ymax) + 
                expand_limits(x=c(0, nrow(sub)+1)) +
                xlab("Event number") + ylab("Normalized intensity") + thms +
                ggtitle(bquote(bold(.(bc_labs[ids == id]))*
                        scriptstyle(" ("*.(n)*" events)")))
        }
        
        # throw informative warning about populations with 
        # no or less than 50 event assignments
        if (!is.null(skipped)) {
            if (length(skipped) == 1) {
                warning("Barcode ID ", paste(skipped), 
                    " has no or less than 50 event assignments;",
                    " no plot has been generated.\n",
                    "It is recommended to remove this population",
                    " or reconsider its separation cutoff.", call.=FALSE)
            } else {
                warning("Barcodes IDs ", paste(skipped, collapse=", "), 
                    " have no or less than 50 event assignments; ",
                    "their plots have not been generated.\n",
                    "It is recommended to remove these populations",
                    " or reconsider their separation cutoffs.", call.=FALSE)
            }
        }
        
        if (length(p) != 0) {
            if (!is.null(out_path)) {
                pdf(file.path(out_path, paste0("event_plot", name_ext, ".pdf")), 
                    width=10, height=5)
                for (i in 1:length(p)) plot(p[[i]])
                dev.off()
            } else {
                for (i in 1:length(p)) plot(p[[i]])
            }
        }
    })