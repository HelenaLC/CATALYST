#' @rdname plotExprs
#' @title Plot expressions
#' 
#' @description 
#' Plots the smoothed densities of arcsinh-transformed marker intensities.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param color_by 
#'   character string. Has to appear as a column name of \code{rowData(x)}.
#'   Specifies the color coding.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return a \code{\link{ggplot}} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotExprs(re)
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
# ------------------------------------------------------------------------------

setMethod(f="plotExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x, color_by="condition") {
        
        df <- data.frame(exprs(x), rowData(x))
        df <- melt(df, id.vars=names(rowData(x)),
            variable.name="antigen", value.name="expression")
        
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
