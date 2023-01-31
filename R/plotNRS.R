#' @rdname plotNRS
#' @title Plot non-redundancy scores
#' 
#' @description Plots non-redundancy scores (NRS) by feature 
#' in decreasing order of average NRS across samples.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param features a character vector specifying 
#'   which antigens to use for clustering; valid values are
#'   \code{"type"/"state"} for \code{type/state_markers(x)} 
#'   if \code{rowData(x)$marker_class} have been specified; 
#'   a subset of \code{rownames(x)}; NULL to use all features.
#' @param color_by character string specifying the color coding;
#'   valid values are \code{namescolData(x))}.
#' @param assay character string specifying which assay data 
#'   to use; valid values are \code{assayNames(x)}.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Crowell HL, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' plotNRS(sce, features = NULL)   # default: all markers
#' plotNRS(sce, features = "type") # type-markers only
#' 
#' @import ggplot2
#' @importFrom Matrix colMeans
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay colData
#' @export

plotNRS <- function(x, features = NULL, 
    color_by = "condition", assay = "exprs") {

    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"), 
        is.character(color_by), length(color_by) == 1,
        sum(color_by == names(colData(x))) == 1)
    
    .check_assay(x, assay)
    features <- .get_features(x, features)
    y <- assay(x, assay)
    y <- y[features, ]
    
    # calculate NRS
    cs_by_s <- split(seq_len(ncol(x)), x$sample_id)
    nrs <- lapply(cs_by_s, function(cs)
        .nrs(y[, cs, drop = FALSE]))
    nrs <- do.call("rbind", nrs)
    
    # warn about excluded markers
    n <- nlevels(x$sample_id)-nrow(nrs)
    if (n > 0) message(
        "Couldn't compute NRS for ", n, " samples\n", 
        "due to an insufficient number of cells.")
    
    # add cell metadata
    m <- match(rownames(nrs), x$sample_id)
    df <- data.frame(nrs, check.names = FALSE)
    df[[color_by]] <- x[[color_by]][m]
    df <- melt(df, 
        id.vars = color_by,
        value.name = "NRS",
        variable.name = "antigen")
    
    # reorder antigens by mean NRS
    avg_nrs <- colMeans(nrs, na.rm = TRUE)
    o <- names(sort(avg_nrs, decreasing = TRUE))
    df$antigen <- factor(df$antigen, levels = o)
    
    ggplot(df, aes_string(x = "antigen", y = "NRS")) +
        geom_point(aes_string(color = color_by), alpha = 0.8,
            position = position_jitter(width = 0.2, height = 0)) +
        geom_boxplot(width = 0.8, fill = NA, outlier.color = NA) +
        stat_summary(fun = "mean", orientation = "x",
            geom = "point", fill = "white", size = 2, shape = 21) +
        guides(col = guide_legend(override.aes = list(alpha = 1, size = 3))) +
        labs(x = NULL, y = "Non-Redundancy score (NRS)") + 
        theme_classic() + theme(
            legend.key.height  =  unit(0.8, "lines"),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_line(color = "grey", linewidth = 0.2))
}