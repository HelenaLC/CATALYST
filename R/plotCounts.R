#' @rdname plotCounts
#' @title Plot cell counts
#' 
#' @description Barplot of the number of cells measured for each sample.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param color_by character string specifying 
#'   a \code{colData(x)} column to color by.
#' @param anno logical. Whether to annotate cell numbers as text above bars.
#' 
#' @author Helena Lucia Crowell
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
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' plotCounts(sce)
#' 
#' @import ggplot2
#' @importFrom methods is
#' @export

plotCounts <- function(x, color_by="condition", anno=TRUE) {
    .check_sce(x)
    .check_cd_factor(x, color_by)
    stopifnot(is.logical(anno), length(anno) == 1)

    p <- ggplot(ei(x), aes_string(x="sample_id", y="n_cells",
        fill=color_by)) + geom_bar(stat="identity", width=.75) +
        coord_cartesian(clip = "off") + theme_minimal() + theme(
            panel.grid.minor=element_blank(), 
            panel.grid.major.x=element_blank(),
            panel.grid.major.y=element_line(color="grey", size=.3), 
            axis.text=element_text(color="black"),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1))
    if (anno) {
        p + geom_label(aes_string(label="n_cells", vjust=-.1), 
            fontface="bold.italic", fill=NA, label.size=0)
    } else {
        p
    }
}
