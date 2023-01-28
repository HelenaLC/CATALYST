#' @rdname plotCounts
#' @title Plot cell counts
#' 
#' @description Barplot of the number of cells measured for each sample.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param group_by character string specifying a non-numeric 
#'   cell metadata column to group by (determines x-axis ticks);
#'   valid values are \code{names(colData(x))}.
#' @param color_by character string specifying a non-numeric
#'   cell metadata column to color by (determines grouping of bars);
#'   valid values are \code{names(colData(x))}; NULL for no color.
#' @param prop logical specifying whether to plot relative abundances 
#'   (frequencies) for each group rather than total cell counts;
#'   bars will be stacked when \code{prop = TRUE} and dodged otherwise.
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
#' 
#' # plot number of cells per sample, colored by condition
#' plotCounts(sce, 
#'   group_by = "sample_id", 
#'   color_by = "condition")
#'   
#' # same as above, but order by patient  
#' plotCounts(sce, 
#'   group_by = "patient_id", 
#'   color_by = "condition")
#' 
#' # total number of cell per patient
#' plotCounts(sce, 
#'   group_by = "patient_id", 
#'   color_by = NULL)
#' 
#' # plot proportion of cells from each patient by condition
#' plotCounts(sce, 
#'   prop = TRUE, 
#'   group_by = "condition", 
#'   color_by = "patient_id")
#' 
#' @import ggplot2
#' @importFrom methods is
#' @export

plotCounts <- function(x, 
    group_by = "condition",
    color_by = group_by, 
    prop = FALSE) {
    
    # check validity of input arguments
    .check_sce(x)
    stopifnot(!is.null(group_by))
    .check_cd_factor(x, color_by)
    .check_cd_factor(x, group_by)
    stopifnot(is.logical(prop), length(prop) == 1)

    df <- data.frame(
        x[[group_by]],
        row.names = NULL, 
        check.names = FALSE)
    if (!is.null(color_by) && color_by != group_by) {
        df[[color_by]] <- x[[color_by]]
        y <- table(df)
        if (prop) 
            y <- y/rowSums(y)
        df <- melt(y, 
            varnames = c(group_by, color_by))
    } else {
        y <- table(df)
        if (prop)
            y <- y/sum(y)
        df <- data.frame(
            value = c(y),
            row.names = NULL)
        df[[group_by]] <- rownames(y)
    }
    
    ggplot(df[df$value != 0, ], aes_string(
        x = group_by, y = "value", fill = color_by)) + 
        geom_bar(stat = "identity", 
            col = ifelse(prop, "white", NA), 
            position = ifelse(prop, "stack", "dodge2")) +
        scale_y_continuous(expand = c(0, 0)) +
        ylab(ifelse(prop, "frequency", "n_cells")) +
        coord_cartesian(clip = "off") + 
        theme_minimal() + theme(
            panel.grid.minor = element_blank(), 
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey", linewidth = 0.2), 
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}
