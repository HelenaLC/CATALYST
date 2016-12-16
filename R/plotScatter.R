# ================================================================================
# Channel-wise medians of measurement vs. compensated data 
# --------------------------------------------------------------------------------

#' @rdname plotScatter
#' @title Channel-wise medians of measurement vs. compensated data
#' 
#' @description For each barcode populations, shows the channel-wise medians
#' computed from measurement and compensated data.
#'
#' @param x        
#' a \code{dbFrame}.
#' @param CM       
#' matrix returned from \code{computeCompmat}.
#' @param out_path
#' a character string. If specified, outputs will be generated in this location. 
#' Defaults to NULL.
#' @param name_ext a character string. If specified, will be appended to the plot's name. Defaults to NULL.
#' 
#' @examples
#' data(ss_exp)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' compMat <- computeCompmat(x = re)
#' plotScatter(x = re, CM = compMat)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette

# --------------------------------------------------------------------------------

plotScatter <- function(x, CM, out_path=NULL, name_ext=NULL) {
    
    ids <- as.numeric(rownames(x@bc_key))
    n_bcs <- length(ids)
    
    nms <- paste(colnames(CM))
    ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
    bc_cols <- which(ms %in% ids)
    bc_range <- min(bc_cols) : max(bc_cols)
    n <- length(bc_range)
    ms <- ms[bc_range]
    
    es <- list(x@exprs, x@exprs %*% CM)
    es <- lapply(es, function(x) x[, bc_range])

    df <- NULL
    for (i in 1:length(es)) {
        for (id in ids) {
            ex <- ms == id
            df <- rbind(df, data.frame(bc=rep(which(ids == id), n-1), 
                m=apply(es[[i]][x@bc_ids == id, !ex], 2, median),
                ch=bc_range[!ex], panel=i))
        }
    }
    df$panel <- factor(df$panel, labels=c("Uncompensated", "Compensated"))
    
    pal <- RColorBrewer::brewer.pal(11, "Spectral")
    if (n > 11) {
        cols <- colorRampPalette(pal)(n)
    } else {
        cols <- sample(pal, n)
    }
    
    p <- ggplot() + theme_bw() + facet_grid(.~panel) +
        geom_point(data=df, aes_string(x="bc", y="m", col="as.factor(ch)"), size=2.5) + 
        scale_colour_manual(values=cols, labels=nms[bc_range], name=NULL) +
        theme(legend.key=element_blank()) + expand_limits(x=c(0, n_bcs+1)) +
        scale_x_discrete(limits=1:n_bcs, labels=nms[bc_cols]) +
        xlab("Barcode population") + ylab("Median counts") +
        theme(panel.border=element_blank(),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_line(size=.2, color="lightgrey"),
            axis.text.x=element_text(vjust=.5, hjust=1, size=5, angle=90),
            axis.text.y=element_text(vjust=.5, hjust=1, size=5))
    
    if (!is.null(out_path)) {
        ggsave(file.path(out_path, paste0("comp_scatter", name_ext,".pdf")), 
            plot=p, width=10, height=5)
    } else {
        plot(p)
    }
}

  