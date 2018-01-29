# ==============================================================================
# Smoothed densities of marker expression
# ------------------------------------------------------------------------------

#' @rdname plotExprs
#' @title Smoothed densities of marker expression
#' 
#' @description 
#' Plots the smoothed densities of arcsinh-transformed marker intensities.
#'
#' @param x a \code{\link{daFrame}}.
#' @param color_by a character string that appears as a column name in the
#' metadata-table of the input \code{daFrame}. Specifies the color coding.
#' In the case of multiple conditions, \code{"condition"} will color points
#' according to each sample's combined conditions.
#' 
#' @return a \code{\link{ggplot}} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotExprs(re)
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
# ==============================================================================

setMethod(f="plotExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x, color_by="condition") {
        df <- data.frame(exprs(x), rowData(x))
        md <- metadata(x)[[1]]
        conds <- grep("condition", colnames(md), value=TRUE)
        if (length(conds) > 1) {
            conds_combined <- apply(md[, conds], 1, paste, collapse="/")
            df$condition <- rep(conds_combined, metadata(x)$n_events)
        }
        df <- melt(df, variable.name="antigen", value.name="expression",
            id.var=c(setdiff(colnames(rowData(x)), "condition"), "condition"))
        ggplot(df, aes_string(x="expression", 
            col=color_by, group="sample_id"), fill=NULL) + 
            facet_wrap(~antigen, ncol=5, scales="free") +
            geom_density() + theme_classic() + theme(
                panel.grid=element_blank(), 
                strip.background=element_blank(),
                strip.text=element_text(face="bold"),
                axis.text=element_text(color="black"), 
                axis.title=element_text(color="black"))
    }
)
