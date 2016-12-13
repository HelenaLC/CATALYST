# ==============================================================================
# Plot spillover matrix heat map
# ------------------------------------------------------------------------------

#' @rdname plotSpillmat
#' @title Spillover matrix heat map
#' 
#' @description 
#' Generates a heat map of the spillover matrix annotated with 
#' estimated spill percentages.
#'
#' @param bc_ms 
#' a vector of numeric masses corresponding to barcode channels.
#' @param CM 
#' matrix returned from \code{computeCompmat}.
#' @param out_path 
#' a character string. If specified, outputs will be generated in this location. 
#' Defaults to NULL.
#' @param name_ext 
#' a character string. If specified, will be appended to the plot's name. 
#' Defaults to NULL.
#' 
#' @examples
#' data(ss_beads)
#' bc_ms <- c(139, 141:156, 158:176)
#' re <- assignPrelim(x = ss_beads, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' compMat <- computeCompmat(x = re)
#' plotSpillmat(bc_ms = bc_ms, CM = compMat)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom grDevices colorRampPalette
#' @export

# ------------------------------------------------------------------------------

plotSpillmat <- function(bc_ms, CM, out_path=NULL, name_ext=NULL) {
    
    nms <- colnames(CM)
    ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
    bc_cols <- which(ms %in% bc_ms)
    inds <- bc_cols - min(bc_cols) + 1
    bc_range <- min(bc_cols) : max(bc_cols)
    n <- length(bc_range)
    SM <- solve(CM)[bc_range, bc_range]
    axis_labs <- paste(colnames(CM)[bc_range])
    axis_labs_cols <- rep("grey", n)
    axis_labs_cols[inds] <- "black"
    
    df <- data.frame(c1=rep(1:n, n), 
        c2=rev(rep(1:n, each=n)), 
        spill=round(100*c(t(SM)), 1))
    
    pal <- colorRampPalette(c("white","lightcoral", "red2", "darkred"))(501)
    spill_labs <- sprintf("%.1f", df$spill)
    spill_labs[df$spill == 0 | df$spill == 100] <- ""
    
    row_labs <- format(round(rowSums(t(matrix(df$spill, nrow=n))) - 100, 2), digits=2)
    col_labs <- format(round(colSums(t(matrix(df$spill, nrow=n))) - 100, 2), digits=2)
    row_labs[-inds] <- col_labs[-inds] <- ""
    
    p <- ggplot(df, aes_string(x="c1", y="c2", fill="spill")) + 
        geom_tile(col="lightgrey") + 
        scale_fill_gradientn(colours=pal, guide=FALSE, limits=c(0, 5), na.value="lightgrey") +
        geom_text(aes_string(label="spill_labs"), col="black", size=2.5) +
        annotate("text", rep(n+1.15, n), 1:n, label=rev(row_labs), fontface="bold", size=2) +
        annotate("text", 1:n, rep(n+1, n), label=col_labs, fontface="bold", size=2) +
        scale_x_discrete(limits=1:n, expand=c(0,0), labels=axis_labs) +
        scale_y_discrete(limits=1:n, expand=c(0,0), labels=rev(axis_labs)) +
        theme_bw() + coord_fixed() + xlab(NULL) + ylab(NULL) +
        theme(panel.grid.major=element_blank(), panel.border=element_blank(),
            axis.text.x=element_text(vjust=.5, hjust=1, size=8, color=axis_labs_cols, angle=90),
            axis.text.y=element_text(vjust=.5, hjust=1, size=8, color=rev(axis_labs_cols)))
    p <- ggplot_gtable(ggplot_build(p))
    p$layout$clip[p$layout$name == "panel"] <- "off"
    
    if (!is.null(out_path)) {
        ggsave(file.path(out_path, paste0("sm_heat", name_ext,".pdf")), 
            plot=p, width=8.1, height=8)
    } else {
        grid.draw(p)
    }
}