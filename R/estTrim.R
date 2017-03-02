# ==============================================================================
# Estimate separation cutoffs
# ------------------------------------------------------------------------------

#' @rdname estTrim
#' @title Estimation of optimal trim value
#' 
#' @description
#' Estimates a trim value that will minimize the sum over squared
#' popultion- and channel-wise squared medians upon compensation.
#'
#' @param x a \code{\link{dbFrame}}.
#' @param min,max,step 
#' specify sequence of trim values for which compensation should be evaluated.
#' @param out_path 
#' specifies in which location output plot is to be generated. Defaults to NULL.
#' @param name_ext 
#' a character string. If specified, will be appended to the output plot's name. 
#' Defaults to NULL.
#' 
#' @return
#' For each value along \code{seq(min, max, step), \code{\link{computeSpillmat}} 
#' will be called with \code{method = "mean"} and the respective trim parameter. 
#' Returned will be the value that minimizes the sum over squared population-
#' wise median counts across all barcodes after compensation.}
#' 
#' @examples
#' data(ss_exp)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' estTrim(x = re, min = 0.06, max = 0.14, step = 0.02)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2
#' @importFrom graphics plot
#' @export

# ------------------------------------------------------------------------------

setMethod(f="estTrim", 
          signature=signature(x="dbFrame"), 
          definition=function(x, min = 0.05, max = 0.20, step = 0.01, 
              out_path = NULL, name_ext = NULL) {
              
              ids <- as.numeric(rownames(x@bc_key))
              nms <- paste(colnames(x@exprs))
              ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
              bc_cols <- which(ms %in% ids)
              bc_range <- min(bc_cols) : max(bc_cols)
              x@exprs <- x@exprs[, bc_range]
              ms <- ms[bc_range]
              
              trim_vals <- seq(min, max, step)
              df <- rss <- NULL
              for (trim in trim_vals) {
                  sm <- computeSpillmat(x=x, method="mean", trim=trim)
                  sm <- make_symetric(sm)
                  tmp <- x@exprs %*% solve(sm)
                  for (id in ids) {
                      median <- apply(tmp[x@bc_ids == id, ms != id], 2, median)
                      df <- rbind(df, data.frame(m=median, trim=trim))
                  }
                  rss <- rbind(rss, data.frame(
                      m=mean((df$m[df$trim == trim])^2), 
                      trim=trim))
              }
              opt <- rss$trim[which.min(rss$m)]
              
              p <- ggplot() + 
                  geom_vline(aes_string(xintercept=opt), lty=3) +
                  geom_jitter(data=df, aes_string(x="trim", y="m"), col="mediumblue", 
                      height=0, width=step/5, size=4, stroke=.5, alpha=.2) + 
                  geom_point(data=rss, aes_string(x="trim", y="m"), size=12, shape=16) +
                  annotate("text", x=rss$trim, y=rss$m, col="white", size=2, 
                      label=paste(round(rss$m, 4))) +
                  scale_x_continuous(breaks=trim_vals, labels=format(trim_vals, 2)) + 
                  xlab("Trim value") + ylab("Median counts") + theme_bw() +
                  theme(legend.key=element_blank(), 
                        panel.border=element_blank(),
                        panel.grid.minor=element_blank(),
                        panel.grid.major=element_line(size=.2, color="lightgrey"),
                        axis.text.x=element_text(vjust=.5, hjust=1))
              
              if (!is.null(out_path)) {
                  ggsave(file.path(out_path, paste0("trim_scatter", name_ext, ".pdf")), plot=p)
              } else {
                  p
              }
              return(opt)
          })
