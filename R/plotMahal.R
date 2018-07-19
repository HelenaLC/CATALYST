#' @rdname plotMahal
#' @title Biaxial plot
#' 
#' @description 
#' Histogram of counts and plot of yields as a function of separation cutoffs.
#'
#' @param x 
#'   a \code{\link{dbFrame}}.
#' @param which 
#'   character string. Specifies which barcode to plot.
#' @param cofactor
#'   numeric. Cofactor used for asinh transformation.
#' @param out_path 
#'   character string. If specified, outputs will be generated here.
#' @param name_ext 
#'   character string. If specified, will be appended to file name. 
#' 
#' @return Plots all inter-barcode interactions for the population specified
#' by argument \code{which}. Events are colored by their Mahalanobis distance. 
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' plotMahal(x = re, which = "B3")
#'
#' @import ggplot2 grid gridExtra
#' @importFrom RColorBrewer brewer.pal
# ------------------------------------------------------------------------------

setMethod(f="plotMahal", 
    signature=signature(x="dbFrame"), 
    definition=function(x, which, cofactor=50, out_path=NULL, name_ext=NULL) {
        
        if (ncol(bc_key(x)) > 7)
            stop("\nToo many barcodes: ",
                "Plotting all inter-barcode interactions is infeasible.\n",
                "Using the default 'mhl_cutoff' value of 30 is recommended.")
        
        inds <- which(bc_ids(x) == which)
        if (length(inds) > 5e3) 
            inds <- sample(inds, 5e3)
        nms <- colnames(exprs(x))
        ms <- as.numeric(get_ms_from_chs(nms))
        bc_cols <- ms %in% colnames(bc_key(x))
        es <- asinh(exprs(x)[inds, bc_cols] / cofactor)
        
        thms <- theme_classic() + theme(
            plot.margin=unit(c(0, 0, 0, 0), "null"),
            axis.title=element_text(size=10), 
            axis.ticks=element_blank(),
            axis.text=element_blank(),
            aspect.ratio=1)
        
        max <- ceiling(max(mhl_dists(x)[inds])/5)*5
        axes_min <- floor(  min(es, na.rm=TRUE)/.1)*.1
        axes_max <- ceiling(max(es, na.rm=TRUE)/.1)*.1
        lims <- c(axes_min, axes_max)
        
        n <- ncol(bc_key(x))
        nPlots <- sum(seq_len(n))
        ps <- vector("list", nPlots)
        
        # histogram of barcode distributions
        hist_inds <- c(1, cumsum(seq(n, 2))+1)
        for (p in hist_inds) {
            i <- j <- match(p, hist_inds)
            df <- data.frame(x=es[, i], y=es[, j], col=mhl_dists(x)[inds])
            ps[[p]] <- ggplot(df) + 
                scale_x_continuous(limits=lims) +
                geom_histogram(aes_string(x="x"), 
                    breaks=seq(axes_min+.05, axes_max-.05, .1),
                    binwidth=.1, fill="black", color=NA) + 
                thms + labs(x=" ", y=" ") + coord_fixed(1)
        }
        
        m <- matrix(NA, n, n)
        m[lower.tri(m, diag=TRUE)] <- seq_along(ps)
        first <- TRUE
        
        # bead vs. bead scatters color coded 
        # according to Mahalanobis distances
        for (i in seq_len(n-1)) {
            for (j in (i+1):n) {
                df <- data.frame(x=es[, i], y=es[, j], col=mhl_dists(x)[inds])
                p <- m[j, i]
                ps[[p]] <- ggplot(df) + labs(x=" ", y=" ") +
                    geom_point(aes_string(x="x", y="y", col="col"), 
                        size=1) + thms + lims(x=lims, y=lims) +
                    guides(colour=FALSE) + scale_color_gradientn(
                        colours=rev(brewer.pal(11, "RdYlBu")),
                        limits=c(0, max), breaks=seq(0, max, 5),
                        name=paste0(which, ": ", 
                            paste(bc_key(x)[which,], collapse="")))
                if (first) {
                    # get color bar
                    lgd <- get_legend(ps[[p]] + 
                            guides(colour=guide_colourbar(
                                title.position="top", title.hjust=.5)) + 
                            theme(legend.direction="horizontal",
                                legend.title=element_text(face="bold"),
                                legend.key.height=unit(1, "line"),
                                legend.key.width=unit(4, "line"),
                                legend.text=element_text(size=8), 
                                legend.key=element_blank()))
                    first <- FALSE
                }
            }
        }
        # add axes labels to left column and bottom row plots
        labs <- colnames(es)
        for (i in seq_len(n)) {
            l <- m[i, 1]
            b <- m[n, i]
            ps[[b]] <- ps[[b]] + labs(x=labs[i])
            ps[[l]] <- ps[[l]] + labs(y=labs[i])
        }
        ps[[nPlots+1]] <- lgd
        m <- rbind(rep(nPlots+1, n), m)
        
        heights <- c(2, rep(5, n))
        widths  <- rep(5, n)
        
        if (!is.null(out_path)) {
            heights <- c(1.5, rep(5, n))
            widths  <- rep(5, n)
            pdf(file.path(out_path, paste0("mahal_plot", name_ext, ".pdf")), 
                width=12, height=12*(sum(heights) / (sum(widths))))
            grid.arrange(grobs=ps, layout_matrix=m, 
                heights=heights, widths=widths)
            dev.off()
        } else {
            heights <- c(2.5, rep(5, n))
            widths  <- rep(5, n)
            grid.arrange(grobs=ps, layout_matrix=m, 
                heights=heights, widths=widths)
        }
    })