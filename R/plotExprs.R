#' @rdname plotExprs
#' @title Expression densities
#' 
#' @description  
#' Plots smoothed densities of marker intensities, with a density curve for 
#' each sample ID, and curves colored by a cell metadata variable of interest.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param features character vector specifying 
#'   which features to invlude; valid values are 
#'   \code{"type"/"state"} for \code{type/state_markers(x)} 
#'   if \code{rowData(x)$marker_class} have been specified; 
#'   a subset of \code{rownames(x)}; NULL to use all features.
#' @param color_by character string specifying 
#'   a non-numeric cell metadata column by which 
#'   to color density curves for each sample; 
#'   valid values are \code{names(colData(x))}.
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

plotExprs <- function(x, features = NULL, color_by = "condition") {
    # check validity of input arguments
    .check_sce(x)
    .check_cd_factor(x, color_by)
    
    # subset features to use
    features <- .get_features(x, features)
    y <- assay(x, "exprs")[features, ]
    
    # construct data.frame include cell metadata
    df <- data.frame(t(y), colData(x), check.names = FALSE)
    gg_df <- melt(df, 
        value.name = "expression",
        variable.name = "antigen", 
        id.vars = names(colData(x)))
    
    ggplot(gg_df, fill = NULL, 
        aes_string(
            x = "expression", y = "..ndensity..",
            col = color_by, group = "sample_id")) + 
        facet_wrap(~ antigen, scales = "free_x") +
        geom_density() + 
        ylab("normalized density") +
        theme_classic() + theme(
            panel.grid = element_blank(), 
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            axis.text = element_text(color = "black"), 
            axis.title = element_text(color = "black"))
}