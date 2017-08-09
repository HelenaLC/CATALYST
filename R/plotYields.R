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
#' @importFrom reshape2 melt
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom htmltools save_html

# ------------------------------------------------------------------------------

setMethod(f="plotYields", 
    signature=signature(x="dbFrame"), 
    definition=function(x, which=0, out_path=NULL, name_ext=NULL) {
        
        ids <- rownames(bc_key(x))
        which <- check_validity_which(which, ids, "yields")
        n_bcs <- nrow(bc_key(x))
        seps <- seq(0, 1, .01)
        bc_labs <- get_bc_labs(x)

        ps <- vector("list", length(which))
        for (i in seq_along(which)) {
            id <- which[i]
            ps[[i]] <- plot_yields(
                id, x, seps, n_bcs, bc_labs)
            if (length(sep_cutoffs(x)) != 0) {
                if (id == 0) {
                    p <- paste0(sprintf("%2.2f", sum(yields(x)[cbind(
                        seq_len(n_bcs), findInterval(sep_cutoffs(x), seps))])
                        /n_bcs*100), "%")
                } else {
                    p <- paste0(sprintf("%2.2f", yields(x)
                        [id, findInterval(sep_cutoffs(x)[id], seps)]*100), "%")
                    ps[[i]] <- ps[[i]] + ggtitle(paste(bc_labs[ids == id], 
                        paste("(separation cutoff ", sep_cutoffs(x)[id],
                        " with", p, "yield)"), sep="<br>"))
                }
            } else if (id != 0) {
                ps[[i]] <- ps[[i]] + ggtitle(bc_labs[ids == id])
            }
            if (id == 0) {
                ps[[i]] <- ggplotly(ps[[i]], tooltip="text") %>% 
                    config(displayModeBar=FALSE)
                ps[[i]] <- hide_legend(ps[[i]])
            } else {
                ps[[i]] <- ggplotly(ps[[i]], 
                    tooltip=c("Cutoff","Yield","Count")) %>% 
                    config(displayModeBar=FALSE)
                ps[[i]] <- layout(ps[[i]], margin=list(t=75))
            }
        }
        if (!is.null(out_path)) {
            outNm <- paste0(paste0("yield_plot", name_ext, ".html"))
            htmltools::save_html(html=ps, file=outNm)
        } else {
            ps
        }
    })
