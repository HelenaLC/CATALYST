#' @rdname plotExprs
#' @title Plot expressions
#' 
#' @description 
#' Plots the smoothed densities of arcsinh-transformed marker intensities.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param color_by character string corresponding to a 
#'   \code{colData(x)} column. Specifies the color coding.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Crowell HL, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return a \code{\link{ggplot}} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' plotExprs(sce)
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment assay colData
#' @export

plotExprs <- function(x, color_by="condition") {
    .check_sce(x)
    stopifnot(is.character(color_by), length(color_by) == 1,
        color_by %in% names(colData(x)))
    
    df <- data.frame(t(assay(x, "exprs")), colData(x))
    df <- melt(df, id.vars=names(colData(x)),
        variable.name="antigen", value.name="expression")
    
    ggplot(df, aes_string(x="expression", 
        col=color_by, group="sample_id"), fill=NULL) + 
        facet_wrap(~antigen, scales="free") +
        geom_density() + theme_classic() + theme(
            panel.grid=element_blank(), 
            strip.background=element_blank(),
            strip.text=element_text(face="bold"),
            axis.text=element_text(color="black"), 
            axis.title=element_text(color="black"))
}