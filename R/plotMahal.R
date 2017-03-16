# ==============================================================================
# Plot counts and yields
# ------------------------------------------------------------------------------

#' @rdname plotMahal
#' @title Biaxial plot
#' 
#' @description 
#' Histogram of counts and plot of yields as a function of separation cutoffs.
#'
#' @param x a \code{\link{dbFrame}}.
#' @param which 
#' specifies which barcode to plot.
#' @param cofactor
#' cofactor used for asinh transformation.
#' @param out_path 
#' a character string. If specified, outputs will be generated in this location. 
#' Defaults to NULL.
#' @param name_ext 
#' a character string. If specified, will be appended to the plot's name. 
#' Defaults to NULL.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' plotMahal(x = re, which = "B3")
#'
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
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
        
        inds <- x@bc_ids == which
        if (sum(inds) > 5e3) 
            inds <- inds[sample(which(inds), 5e3)]
        ms <- gsub("[[:alpha:][:punct:]]", "", colnames(x@exprs))
        es <- asinh(x@exprs[inds, ms %in% colnames(x@bc_key)] / cofactor)
        
        thms <- theme_classic() + theme(
            plot.margin=unit(c(0,0,0,0), "null"),
            axis.title=element_text(size=10), 
            axis.ticks = element_blank(),
            axis.text=element_blank(),
            aspect.ratio=1)
        
        max_dist <- ceiling(max(x@mhl_dists[inds])/5)*5
        max_int  <- max(es)
        
        n <- ncol(x@bc_key)
        p <- list()
        for (i in 1:n) {
            for (j in i:n) {
                df <- data.frame(x=es[, i], y=es[, j], col=x@mhl_dists[inds])
                if (i == j) {
                    p[[length(p) + 1]] <- ggplot(df) + 
                        scale_x_continuous(limits=c(0, max_int)) +
                        geom_histogram(aes_string(x="x"), bins=100,
                            fill="black", color=NA) + 
                        thms + labs(x=" ", y=" ") + coord_fixed(1)
                } else {
                    p[[length(p) + 1]] <- ggplot(df) + 
                        labs(x=" ", y=" ") + thms +
                        scale_x_continuous(limits=c(0, max_int)) +
                        scale_y_continuous(limits=c(0, max_int)) +
                        geom_point(aes_string(x="x", y="y", col="col"), 
                            size=1) +
                        guides(colour=guide_colourbar(
                            title.position="top", 
                            title.hjust=.5)) +
                        scale_color_gradientn(
                            colours=rev(brewer.pal(11, "RdYlBu")),
                            limits=c(0, max_dist), breaks=seq(0, max_dist, 5),
                            name=paste0(which, ": ", 
                                paste(x@bc_key[which,], collapse="")))
                }
                if (i == 1 & j == 2) {
                    lgd <- get_legend(p[[length(p)]] + theme(
                        legend.direction="horizontal",
                        legend.title=element_text(face="bold"),
                        legend.key.height=unit(1, "line"),
                        legend.key.width=unit(4, "line"),
                        legend.text=element_text(size=8), 
                        legend.key=element_blank()))
                }
                p[[length(p)]] <- p[[length(p)]] + guides(colour=FALSE)
                if (i == 1) {
                    p[[length(p)]] <- p[[length(p)]] + ylab(colnames(es)[j])
                }
                if (j == n) {
                    p[[length(p)]] <- p[[length(p)]] + xlab(colnames(es)[i])
                }
            }
        }
        
        m <- matrix(NA, n, n)
        m[lower.tri(m, diag=TRUE)] <- 1:length(p)
        p[[length(p) + 1]] <- lgd
        m <- rbind(rep(length(p), n), m)
        
        heights <- c(2, rep(5, n))
        widths  <- rep(5, n)
        
        if (!is.null(out_path)) {
            heights <- c(1.5, rep(5, n))
            widths  <- rep(5, n)
            pdf(file.path(out_path, paste0("mahal_plot", name_ext, ".pdf")), 
                width=12, height=12*(sum(heights) / (sum(widths))))
            grid.arrange(grobs=p, layout_matrix=m, 
                heights=heights, widths=widths)
            dev.off()
        } else {
            heights <- c(2.5, rep(5, n))
            widths  <- rep(5, n)
            grid.arrange(grobs=p, layout_matrix=m, 
                heights=heights, widths=widths)
        }
    })

