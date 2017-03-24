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
#' @param SM 
#' spillover matrix returned from \code{computeSpillmat}.
#' @param annotate
#' logical. If TRUE (default), spill percentages are shown inside bins and
#' rows/columns are annotated with the total amount of spill caused/received. 
#' @param palette
#' an optional vector of colors to interpolate.
#' 
#' @return plots estimated spill percentages as a heat map. 
#' Colours are ramped to the highest spillover value present
#' 
#' @examples
#' # get single-stained control samples
#' data(ss_exp)
#' 
#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#' 
#' re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' spillMat <- computeSpillmat(x = re)
#' plotSpillmat(bc_ms = bc_ms, SM = spillMat)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom grDevices colorRampPalette
#' @export

# ==============================================================================

plotSpillmat <- function(bc_ms, SM, annotate=TRUE, palette=NULL) {
    
    nms <- colnames(SM)
    ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
    bc_cols <- which(ms %in% bc_ms)
    bc_range <- min(bc_cols) : max(bc_cols)
    SM <- make_symetric(SM)[bc_range, bc_range]
    n <- length(bc_range)
    axis_labs <- nms[bc_range]
    lab_cols <- rep("grey", n)
    lab_cols[axis_labs %in% nms[bc_cols]] <- "black"
    
    df <- data.frame(c1=rep(1:n, n), 
        c2=rev(rep(1:n, each=n)), 
        spill=round(100*c(t(SM)), 1))
    
    max <- ceiling(max(df$spill[df$spill != 100])/.25)*.25
    
    if (is.null(palette)) 
        palette <- c("white","lightcoral", "red2", "darkred")
    if (!palette[1] %in% c("white", "#FFFFFF"))
        palette <- c("white", palette)
    pal <- colorRampPalette(palette)(max*100)
    
    p <- ggplot(df, aes_string(x="c1", y="c2")) + 
        geom_tile(aes_string(fill="spill"), col="lightgrey") + 
        scale_fill_gradientn(colors=pal, limits=c(0, max), guide=FALSE) +
        scale_x_discrete(limits=1:n, expand=c(0,0), labels=axis_labs) +
        scale_y_discrete(limits=1:n, expand=c(0,0), labels=rev(axis_labs)) +
        coord_fixed() + xlab(NULL) + ylab(NULL) + theme_bw() + theme(
            panel.grid.major=element_blank(), panel.border=element_blank(),
            axis.text.x=element_text(vjust=.5, angle=90),
            axis.text.y=element_text(vjust=.5, color=rev(lab_cols)))
    
    if (annotate) {
        spill_labs <- sprintf("%.1f", df$spill)
        spill_labs[df$spill == 0 | df$spill == 100] <- ""
        
        row_sums <- round(rowSums(t(matrix(df$spill, n)))-100, 2)
        col_sums <- round(colSums(t(matrix(df$spill, n)))-100, 2)
        row_labs <- format(row_sums, digits=2)
        col_labs <- format(col_sums, digits=2)
        
        p <- p + geom_text(aes_string(label="spill_labs"), size=3) +
            annotate("text", rep(n+1.15, n), 1:n, label=rev(row_labs), 
                fontface="bold", size=2.5, col=rev(lab_cols)) +
            annotate("text", 1:n, rep(n+1, n), label=col_labs,      
                fontface="bold", size=2.5) 
    }
    
    p <- ggplot_gtable(ggplot_build(p))
    p$layout$clip[p$layout$name == "panel"] <- "off"
    grid::grid.draw(p)
}

