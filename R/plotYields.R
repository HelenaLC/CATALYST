# ==============================================================================
# Plot distribution of barcode separations & 
# yields as a function of separation cutoffs
# ------------------------------------------------------------------------------

#' @rdname plotYields
#' @title Yield plot
#' 
#' @description 
#' Distribution of barcode separations and 
#' yields as a function of separation cutoffs.
#'
#' @param x a \code{\link{dbFrame}}.
#' @param which 0, numeric or character. Specifies which barcode(s) to plot. 
#' Valid values are IDs that occur as row names in the \code{bc_key} of the 
#' supplied \code{\link{dbFrame}}; 0 (the default) will generate a summary plot 
#' with all barcodes.
#' @param annotate
#' logical. If TRUE (default) and the \code{sep_cutoffs} slot of the supplied
#' \code{\link{dbFrame}} is not empty, vertical lines will be drawn at cutoff
#' values and the resulting yield will be included in the plot title.
#' @param legend
#' logical. Specifies if a legend should be included. 
#' This will only affect the summary plot (\code{which=0}).
#' @param out_path 
#' a character string. If specified, outputs will be generated 
#' in this location. Defaults to NULL.
#' @param name_ext 
#' a character string. If specified, will be appended to the plot's name. 
#' Defaults to NULL.
#'
#' @details
#' The overall yield that will be achieved upon application of the specified 
#' set of separation cutoffs is indicated in the summary plot. Respective 
#' separation thresholds and their resulting yields are included in each 
#' barcode's plot. The separation cutoff value should be chosen such that
#' it appropriately balances confidence in barcode assignment and cell yield.
#' 
#' @return plots the distribution of barcode separations and yields upon 
#' debarcoding as a function of separation cutoffs. If available, currently 
#' used separation cutoffs as well as their resulting yields will be indicated 
#' in the plot`s main title.
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
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom stats predict smooth.spline
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette pdf dev.off

# ------------------------------------------------------------------------------

setMethod(f="plotYields", 
    signature=signature(x="dbFrame"), 
    definition=function(x, which=0, annotate=TRUE, legend=TRUE, 
        out_path=NULL, name_ext=NULL) {
        
        ids <- rownames(bc_key(x))
        which <- check_validity_which(which, ids, "yields")
        n_bcs <- nrow(bc_key(x))
        seps <- seq(0, 1, .01)
        bc_labs <- get_bc_labs(x)

        ps <- list()
        for (id in which) {
            ps[[length(ps)+1]] <- plot_yield(
                id, x, seps, n_bcs, legend, bc_labs)
            if (annotate && length(sep_cutoffs(x)) != 0) {
                if (id == 0) {
                    p <- paste0(sprintf("%2.2f", sum(yields(x)[cbind(1:n_bcs,
                        match(sep_cutoffs(x), seps))])/n_bcs*100), "%")
                    ps[[length(ps)]] <- ps[[length(ps)]] + 
                        ggtitle(bquote(bold(.(p))*" overall yield")) 
                } else {
                    p <- paste0(sprintf("%2.2f", yields(x)
                        [id, seps %in% sep_cutoffs(x)[id]]*100), "%")
                    ps[[length(ps)]] <- ps[[length(ps)]] + ggtitle(
                        bquote(bold(.(bc_labs[ids == id]))*scriptstyle(
                        " (separation cutoff "*.(sep_cutoffs(x)[id])
                        *" with "*.(p)*" yield)")))
                }
            }
        }

        if (!is.null(out_path))
            pdf(file.path(out_path, paste0("yield_plot", name_ext, ".pdf")),
                height=6, width=12)
        suppressWarnings(lapply(ps, plot))
        if (!is.null(out_path))
            dev.off()
    })

