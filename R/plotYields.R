# ==============================================================================
# Plot x@counts and x@yields
# ------------------------------------------------------------------------------

#' @rdname plotYields
#' @title Yield plot
#' 
#' @description 
#' Histogram of x@counts and plot of x@yields as a function of separation cutoffs.
#'
#' @param x a \code{\link{dbFrame}}.
#' @param which 
#' numeric. Specifies which barcode to plot. 
#' 0, numeric or character. Specifies which barcode(s) to plot. 
#' Valid values are IDs that occur as row names in the \code{bc_key} 
#' of the supplied \code{\link{dbFrame}}; 
#' 0 (default) will generate a summary plot with all barcodes.
#' @param annotate
#' logical. If TRUE (default) and the \code{sep_cutoffs} slot of the supplied
#' \code{\link{dbFrame}} is not empty, vertical lines will be drawn at cutoff
#' values and the resulting yield will be included in the plot title.
#' @param legend
#' logical. Specifies if a legend should be included. 
#' This will only affect the summary plot (\code{which=0}).
#' @param out_path 
#' a character string. If specified, outputs will be generated in this location. 
#' Defaults to NULL.
#' @param name_ext 
#' a character string. If specified, will be appended to the plot's name. 
#' Defaults to NULL.
#'
#' @details 
#' The overall yield that will be achieved upon application of the specified 
#' separation cutoffs is indicated in the summary plot. Respective separation 
#' thresholds and their resulting x@yields are included in each barcode's plot. 
#' The separation cutoff value should be chosen such that it appropriately 
#' balances confidence in barcode assignment and cell yield.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' re <- estCutoffs(x = re)
#' # all samples summary plot
#' plotYields(x = re, which = 0)
#' # plot for specific sample
#' plotYields(x = re, which = "D1")
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
              
              which <- check_validity_which(
                  which, rownames(x@bc_key), "x@yields")
              
              ids <- rownames(x@bc_key)
              n_bcs <- length(ids)
              nms <- colnames(x@exprs)
              ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
              bc_chs <- which(ms %in% as.numeric(colnames(x@bc_key)))
              
              # get barcode labels: 
              # channel name if barcodes are single-positive, 
              # barcode ID and binary code elsewise 
              if (sum(rowSums(x@bc_key) == 1) == n_bcs) {
                  bc_labs <- paste(colnames(x@exprs))[bc_chs]
              } else {
                  bc_labs <- paste0(rownames(x@bc_key), ": ", 
                      apply(x@bc_key, 1, function(x) paste(x, collapse="")))
              }
              
              seps <- seq(0, 1, .01)
              n_seps <- length(seps)
              
              # plotting aesthetics
              pal <- RColorBrewer::brewer.pal(11, "Spectral")
              if (n_bcs > 11) {
                  cols <- colorRampPalette(pal)(n_bcs)
              } else {
                  cols <- sample(pal, n_bcs)
              }
              tick_labs <- paste(round(seps[seq(1, n_seps, 10)], 1))
              thms <- theme_bw() + theme(
                  panel.grid.major=element_line(color="lightgrey"),
                  panel.grid.minor=element_blank())
              
              # initialize list of count histograms (h) and yield plots (l)
              h <- l <- list()
              for (i in which) {
                  if (i == 0) {
                      # summary plot with all barcodes
                      yield <- 0
                      yield <- sum(sapply(1:n_bcs, 
                          function(k) yield <- sum(yield, 
                              x@yields[k, seps == x@sep_cutoffs[k]]))) / n_bcs
                      yield <- paste0(round(yield*100, 2), "%")
                      
                      df_h <- data.frame(sep=seps, count=colSums(x@counts))
                      df_l <- data.frame(sep=rep(seps, n_bcs), 
                          yield=c(t(x@yields)), bc=rep(1:n_bcs, each=n_seps))
                      
                      # barcode-specific x@yields plot
                  } else {
                      ind <- ids == i 
                      df_h <- data.frame(sep=seps, count=x@counts[ind, ])
                      df_l <- data.frame(sep=seps, yield=x@yields[ind, ])
                      df_c <- data.frame(sep=seps, pred=predict(
                          smooth.spline(df_l$sep, df_l$yield))$y)
                  }
                  
                  # x@counts as a function of barcode separation
                  h[[length(h)+1]] <- ggplot() + 
                      geom_bar(data=df_h, aes_string(x="sep", y="count"), 
                               stat="identity", width=1/n_seps, 
                               fill="navy", colour="cornflowerblue") +
                      scale_x_continuous(breaks=seq(0, 1, .1)) +
                      scale_y_continuous(labels=scientific_10) +
                      xlab("Barcode separation") + ylab("Event count") + thms
                  
                  # x@yields as a function of separation cutoff
                  l[[length(l)+1]] <- ggplot() + 
                      scale_x_continuous(breaks=seq(0, 1, .1)) +
                      scale_y_continuous(breaks=seq(0, 1, .25), 
                          labels=paste0(seq(0, 1, .25)*100, "%")) +
                      xlab("Barcode separation threshold") + 
                      ylab("Yield after debarcoding") + thms
                  
                  if (i == 0) {
                      # scale colour according to barcode
                      l[[length(l)]] <- l[[length(l)]] +
                          geom_line(data=df_l, size=.5, 
                              aes_string(x="sep", y="yield", col="as.factor(bc)")) +
                          scale_colour_manual(values=cols, name=NULL, labels=bc_labs)
                      if (legend) {
                          l[[length(l)]] <- l[[length(l)]] +
                              guides(colour=guide_legend(override.aes=list(size=1)))
                          # store legend for later and remove
                          lgd <- get_legend(l[[length(l)]] + theme(
                              legend.text=element_text(size=8), 
                              legend.key=element_blank()))
                      } 
                      l[[length(l)]] <- l[[length(l)]] + guides(colour=FALSE)
                  } else {
                      # add vertical line at cutoff estimate
                      h[[length(h)]] <- h[[length(h)]] +
                          geom_vline(aes_string(xintercept="x@sep_cutoffs[ind]"), lty=3)
                      
                      # add vertical line at cutoff estimate and smooth spline fit
                      l[[length(l)]] <- l[[length(l)]] +
                          geom_point(data=df_l, aes_string(x="sep", y="yield"), 
                                     pch=21, size=3, alpha=.5, stroke=.5,
                                     col="navy", fill="cornflowerblue") +
                          geom_line(data=df_c, aes_string(x="sep", y="pred"), col="red") +
                          geom_vline(aes_string(xintercept="x@sep_cutoffs[ind]"), lty=3)
                  }
                  
                  # adjust widths for grid layout
                  gp1 <- ggplot_gtable(ggplot_build(h[[length(h)]]))
                  gp2 <- ggplot_gtable(ggplot_build(l[[length(l)]]))
                  maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
                  gp1$widths[2:3] <- gp2$widths[2:3] <- maxWidth
                  h[[length(h)]] <- gp1
                  l[[length(l)]] <- gp2
              }
              
              if (0 %in% which) {
                  # no title if 'sep_cutoffs' are missing or 'annotate' = FALSE,
                  # overall yield in % elsewise 
                  if (length(x@sep_cutoffs) == 0 | !annotate) {
                      main <- NULL
                  } else {
                      main <- textGrob(bquote(
                          bold(.(yield))*" overall yield"), just="top")
                  }
                  
                  if (!is.null(out_path))
                      pdf(file.path(out_path, 
                          paste0("summary_yield_plot", name_ext, ".pdf")), 
                          width=10, height=6)
                  ind <- which(which == 0)
                  if (legend) {
                      grid.arrange(h[[ind]], l[[ind]], lgd, ncol=2, nrow=2, 
                          layout_matrix=rbind(c(1, 3), c(2, 3)), 
                          widths=c(8, 2), heights=c(4, 4),
                          top=main)
                  } else {
                      grid.arrange(h[[ind]], l[[ind]], nrow=2, 
                          widths=10, heights=c(4, 4),
                          top=main)
                  }
                  if (!is.null(out_path)) dev.off()
              }
              
              if (sum(fil <- which != 0) != 0) {
                  if (!is.null(out_path))
                      pdf(file.path(out_path, 
                          paste0("yield_plot", name_ext, ".pdf")), 
                          width=8, height=6)
                  for (i in which[fil]) {
                      ind <- match(i, rownames(x@bc_key))
                      j <- which(which == i)
                      perc <- paste0(round(x@yields[ind, sapply(seps, function(x) 
                          isTRUE(all.equal(x, x@sep_cutoffs[ind])))], 3)*100, "%")
                      grid.arrange(h[[j]], l[[j]], nrow=2, widths=8, heights=c(3, 3),
                                   top=textGrob(bquote(bold(.(bc_labs[ind]))*scriptstyle(
                                       " (current cutoff "*.(x@sep_cutoffs[ind])*" with "*.(perc)*" yield)")), just="top"))
                  }
                  if (!is.null(out_path)) dev.off()
              }
          })