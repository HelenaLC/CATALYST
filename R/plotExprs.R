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
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2
#' @importFrom flowCore exprs fsApply parameters
# ==============================================================================

setMethod(f="plotExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x) {
        # arcsinh-transformation & column subsetting
        d <- parameters(data(x)[[1]])$desc
        inds <- d %in% panel(x)$Antigen[panel(x)$Lineage | panel(x)$Functional]
        es <- fsApply(fs, function(ff) asinh(exprs(ff)[, inds]/5))
        # set antigens only as column names
        colnames(es) <- d[inds]
        n_events <- fsApply(data(x), nrow)
        df <- data.frame(es,
            sample_id=rep(md$sample_id, n_events),
            condition=rep(md$condition, n_events))
        df <- melt(df, id.var=c("sample_id", "condition"), 
            variable.name="antigen", value.name="expression")
        ggplot(df, aes_string(x="expression", 
            col="condition", group="sample_id")) + 
            facet_wrap(~antigen, ncol=5, scales="free") + 
            geom_density() + theme_classic() + theme(
                panel.grid=element_blank(), 
                strip.background=element_blank(),
                strip.text=element_text(face="bold"),
                axis.text=element_text(color="black"), 
                axis.title=element_text(color="black"))
    }
)