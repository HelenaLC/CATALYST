# ==============================================================================
# Smoothed densities of marker expression
# ------------------------------------------------------------------------------

#' @rdname plotExprs
#' @title Smoothed densities of marker expression
#' 
#' @description 
#' Sample-wise smoothed densities of marker expressions.
#'
#' @param x a \code{\link{daFrame}}.
#' 
#' @return Plots sample-wise smoothed densities of marker expressions.
#'
#' @references 
#' Nowicka M, Krieg C, Weber LM et al.
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotExprs(re)
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2
#' @importFrom reshape2 melt
# ==============================================================================

setMethod(f="plotExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x) {
        df <- data.frame(exprs(x), 
            sample_id=sample_ids(x), 
            condition=conditions(x))
        df <- melt(df, id.var=c("sample_id", "condition"),
            variable.name="antigen", value.name="expression")
        ggplot(df, aes_string(x="expression", 
            col="condition", group="sample_id"), fill=NULL) + 
            facet_wrap(~antigen, ncol=5, scales="free") +
            geom_density() + theme_classic() + theme(
                panel.grid=element_blank(), 
                strip.background=element_blank(),
                strip.text=element_text(face="bold"),
                axis.text=element_text(color="black"), 
                axis.title=element_text(color="black"))
    }
)