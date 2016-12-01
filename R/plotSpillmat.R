# ================================================================================
# Plot spillover matrix heat map
# --------------------------------------------------------------------------------

#' @rdname plotSpillmat
#' @title Spillover matrix heat map
#' 
#' @description Generates a heat map of the spillover matrix annotated with estimated spill percentages.
#'
#' @param bc_ms    a vector of numeric masses corresponding to barcode channels.
#' @param CM       matrix returned from \code{computeCompmat}.
#' @param out_path specifies in which location output plot is to be generated.
#' @param name_ext a character string. If specified, will be appended to the plot's name. Defaults to NULL.
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

# --------------------------------------------------------------------------------

plotSpillmat <- function(bc_ms, CM, out_path=NULL, name_ext=NULL) {
    
    # get barcode columns
    nms <- colnames(CM)
    ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
    bc_cols <- which(ms %in% bc_ms)
    inds <- bc_cols - min(bc_cols) + 1
    bc_range <- min(bc_cols) : max(bc_cols)
    n <- length(bc_range)
    SM <- solve(CM)[bc_range, bc_range]
    diag(SM) <- 1
    
    labs <- paste(colnames(CM)[bc_range])
    labs_cols <- rep("grey", n)
    labs_cols[inds] <- "black"
    
    df <- data.frame(c1=rep(1:n, n), 
        c2=rev(rep(1:n, each=n)), 
        spill=round(100*c(t(SM)), 2))
    
    vals <- unique(df$spill)
    min_val <- floor(min(df$spill))
    max_val <- max(vals[!(vals == 100)])
    n_cols <- diff(c(min_val, max_val)) * 100

    pal <- colorRampPalette(c("white", "lightcoral", "red2", "darkred"))(n_cols)
    
    anno <- df$spill
    anno[df$spill == 0] <- ""
    
    anno_cols <- matrix(rep(rep("black", n), n), ncol=n)
    diag(anno_cols) <- "white"
    anno_cols <- c(t(anno_cols))
    
    p <- ggplot(df, aes_string(x="c1", y="c2", fill="spill")) + 
        geom_tile(col="lightgrey") + 
        scale_fill_gradientn(colours=pal, guide=FALSE, limits=c(min_val, max_val)) +
        geom_text(aes_string(label="anno"), size=2, col=anno_cols) +
        scale_x_discrete(limits=1:n, labels=labs) +
        scale_y_discrete(limits=1:n, labels=rev(labs)) +
        theme_bw() + coord_fixed() + xlab(NULL) + ylab(NULL) +
        theme(panel.grid.major=element_blank(), panel.border=element_blank(),
            axis.text.x=element_text(vjust=.5, hjust=1, size=8, color=labs_cols, angle=90),
            axis.text.y=element_text(vjust=.5, hjust=1, size=8, color=rev(labs_cols)))
    
    if (!is.null(out_path)) {
        ggsave(file.path(out_path, paste0("sm_heat", name_ext,".pdf")), 
            plot=p, width=8, height=8)
    } else {
        p
    }
}