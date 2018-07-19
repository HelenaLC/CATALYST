#' @rdname plotSpillmat
#' @title Spillover matrix heat map
#' 
#' @description Generates a heat map of the spillover matrix 
#' annotated with estimated spill percentages.
#'
#' @param bc_ms 
#'   a vector of numeric masses corresponding to barcode channels.
#' @param SM 
#'   spillover matrix returned from \code{computeSpillmat}.
#' @param out_path 
#'   character string. If specified, outputs will be generated here.
#' @param name_ext 
#'   character string. If specified, will be appended to the plot's name. 
#' @param annotate
#'   logical. If TRUE (default), spill percentages are shown inside bins 
#'   and rows are annotated with the total amount of spill received. 
#' @param plotly
#'   logical. Should an interactive plot be rendered?
#' @param isotope_list
#'   named list. Used to validate the input spillover matrix.
#'   Names should be metals; list elements numeric vectors of their isotopes.
#'   See \code{\link{isotope_list}} for the list of isotopes used by default.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @return Plots estimated spill percentages as a heat map. 
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
#' @import ggplot2 grid gridExtra
#' @importFrom grDevices colorRampPalette
#' @importFrom htmltools save_html
#' @importFrom plotly ggplotly
# ------------------------------------------------------------------------------

setMethod(f="plotSpillmat",
    signature=signature(bc_ms="numeric", SM="matrix"),
    definition=function(bc_ms, SM, out_path=NULL, 
        name_ext=NULL, annotate=TRUE, plotly=TRUE, 
        isotope_list=CATALYST::isotope_list) {
    
        SM <- check_sm(SM, isotope_list)
        nms <- colnames(SM)
        ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
        bc_cols <- which(ms %in% bc_ms)
        bc_range <- min(bc_cols) : max(bc_cols)
        SM <- make_symetric(SM)[bc_range, bc_range]
        n <- length(bc_range)
        axis_labs <- nms[bc_range]
        ex <- !axis_labs %in% nms[bc_cols]
        lab_cols <- rep("black", n)
        lab_cols[ex] <- "grey"
        
        df <- reshape2::melt(100*SM)
        colnames(df) <- c("Emitting", "Receiving", "spill")
        df$Spillover <- paste0(sprintf("%2.3f", df$spill), "%")
        max <- ceiling(max(100*SM[row(SM) != col(SM)])/.25)*.25
        overallReceived <- paste0(sprintf("%2.2f", 
            rowSums(t(matrix(df$spill, n)))-100, "%"))
        overallReceived[ex] <- NA
        
        p <- ggplot(df, 
            aes_string(x="Receiving", y="Emitting", group="Spillover")) + 
            geom_tile(aes_string(fill="spill"), col="lightgrey") + 
            scale_fill_gradientn(
                colors=c("white", "lightcoral", "red2", "darkred"), 
                limits=c(0, max), na.value="lightgrey", guide=FALSE) +
            scale_x_discrete(limits=colnames(SM), expand=c(0,0)) +
            scale_y_discrete(limits=rev(rownames(SM)), expand=c(0,0)) +
            coord_fixed() + labs(x=NULL, y=NULL) + theme_bw() + theme(
                panel.grid=element_blank(), panel.border=element_blank(),
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
        
        if (annotate) {
            spill_labs <- sprintf("%.1f", df$spill)
            spill_labs[df$spill == 0 | df$spill == 100] <- ""
            if (is.null(out_path)) size <- 3 else size <- 2
            p <- p + geom_text(aes_string(label="spill_labs"), size=size)
        }
        
        if (plotly)
            p <- ggplotly(p, width=720, height=720,
                tooltip=c("Emitting", "Receiving", "Spillover")) %>%
                layout(margin=list(l=72, b=58)) 
        
        if (!is.null(out_path)) {
            if (class(p)[1] == "plotly") {
                htmltools::save_html(p, 
                    file.path(out_path, paste0("SpillMat", name_ext, ".html")))
            } else {
                ggsave(plot=p, width=7, height=7,
                    file.path(out_path, paste0("SpillMat", name_ext, ".pdf")))
            }
        } else {
            p
        }
    }
)

# ------------------------------------------------------------------------------
#' @rdname plotSpillmat
setMethod(f="plotSpillmat",
    signature=signature(bc_ms="ANY", SM="data.frame"),
    definition=function(bc_ms, SM, 
        out_path=NULL, name_ext=NULL, annotate=TRUE, plotly=TRUE) {
        plotSpillmat(bc_ms, as.matrix(SM), out_path=NULL, 
            name_ext=NULL, annotate=TRUE, plotly=TRUE)
    }
)

# ------------------------------------------------------------------------------
#' @rdname plotSpillmat
setMethod(f="plotSpillmat",
    signature=signature(bc_ms="character", SM="ANY"),
    definition=function(bc_ms, SM, 
        out_path=NULL, name_ext=NULL, annotate=TRUE, plotly=TRUE) {
        plotSpillmat(as.numeric(bc_ms), SM, out_path=NULL, 
            name_ext=NULL, annotate=TRUE, plotly=TRUE)
    }
)