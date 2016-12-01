# ================================================================================
# Plot events
# --------------------------------------------------------------------------------

#' @rdname plotEvents
#' @title Event plot
#' 
#' @description 
#' Shows normalized barcode intensities for a given barcode.
#'
#' @param x        a \code{\link{dbFrame}} ccontaining normalized barcode intensities, 
#'                 the debarcoding key, and barcode IDs.
#' @param which_bc numeric or "all". Specifies which barcode to plot. Defaults to "all".
#' @param n_events numeric or "all". Specifies how many events to plot per barcode. Defaults to 100.
#' @param out_path character string. Specifies in which location output plot is to be generated.
#' @param name_ext a character string. If specified, will be appended to the plot's name. Defaults to NULL.
#' 
#' @details 
#' EXPLAIN NORMALIZATION HERE
#' 
#' @examples
#' data(ss_beads)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_beads, y = bc_ms)
#' plotEvents(x = re, which_bc = 139)
#' 
#' @references
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette pdf dev.off

# --------------------------------------------------------------------------------

setMethod(f="plotEvents", 
    signature=signature(x="dbFrame"), 
    definition=function(x, which_bc="all", n_events=100, out_path=NULL, name_ext=NULL) {
        
        normed_bcs <- x@normed_bcs
        bc_key <- x@bc_key
        bc_ids <- x@bc_ids
        nms <- colnames(x@exprs)
        ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
        
        ids <- sort(unique(bc_ids))
        n_ids <- length(which(ids!=0))
        n_bcs <- nrow(bc_key)
        
        if (sum(rowSums(bc_key) == 1) == n_bcs) {
            bc_labs <- paste(colnames(x@exprs))[
                ms %in% as.numeric(colnames(bc_key))]
        } else {
            bc_labs <- paste(apply(bc_key, 1, function(x) paste(x, collapse="")))
        }
        if (0 %in% ids) 
            bc_labs <- c("Unassigned", bc_labs)
        
        labs <- colnames(normed_bcs)
        pal <- RColorBrewer::brewer.pal(11,"Spectral")
        if (n_ids > 11) {
            cols <- colorRampPalette(pal)(ncol(normed_bcs))
        } else {
            cols <- sample(pal,ncol(normed_bcs))
        }
        aes <- theme_bw() + theme(legend.key=element_blank(),
            panel.grid.major=element_line(color="lightgrey"),
            panel.grid.minor=element_blank())
        
        if ("all" %in% which_bc) {
            which_ids <- ids
        } else {
            which_ids <- which_bc
        }
        
        p <- list()
        for (id in which_ids) {
            inds <- which(bc_ids == id)
            n <- length(inds)
            if (length(inds) == 0) next
            if (is.numeric(n_events) & length(inds) > n_events)
                inds <- sort(sample(inds, n_events))
            psize <- 1 + 100 / length(inds)
            sub <- normed_bcs[inds, ]
            
            tcks <- c(1,nrow(sub))
            ymin <-   floor(4 * min(sub)) / 4
            ymax <- ceiling(4 * max(sub)) / 4
            
            df <- data.frame(event=rep(1:nrow(sub), each=ncol(sub)),
                intensity=c(t(sub)),
                bc=rep(1:ncol(sub), nrow(sub)))

            p[[length(p) + 1]] <- ggplot(df) + geom_point(stroke=.5, size=psize, 
                aes_string(x="event", y="intensity", col="as.factor(bc)", alpha=.8)) +
                scale_colour_manual(values=cols, name=NULL, labels=labs) + 
                guides(alpha=FALSE, colour=guide_legend(override.aes=list(size=3))) +
                scale_x_discrete(limits=NULL, labels=NULL) + ylim(ymin, ymax) + 
                expand_limits(x=c(0, nrow(sub)+1)) +
                xlab("Event number") + ylab("Normalized intensity") + aes +
                ggtitle(bquote(bold(.(bc_labs[ids == id]))*
                        scriptstyle(" ("*.(n)*" events)")))
        }
        
        if (!is.null(out_path)) {
            pdf(file.path(out_path, paste0("ep", name_ext, ".pdf")), 
                width=10, height=5)
            for (i in 1:length(p)) plot(p[[i]])
            dev.off()
        } else {
            for (i in 1:length(p)) plot(p[[i]])
        }
    })