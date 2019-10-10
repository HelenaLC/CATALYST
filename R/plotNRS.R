#' @rdname plotNRS
#' @title Plot non-redundancy scores
#' 
#' @description 
#' Plots non-redundancy scores (NRS) by features in decreasing order.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param features
#'   character string specifying which features to include. Defaults to NULL 
#'   (= all features). Alternatively, if the \code{rowData(x)$marker_class} 
#'   column is specified, can be one of "type", "state", or "none".
#' @param color_by 
#'   character string. Has to appeara as a column name of \code{colData(x)}.
#'   Specifies the color coding.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
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
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay
#' @export
# ------------------------------------------------------------------------------

plotNRS <- function(x, features=NULL, color_by="condition") {
    
    stopifnot(is.character(color_by), 
        color_by %in% colnames(colData(x)))
    
    if (is.null(features)) {
        features <- rownames(x)
    } else if (length(features) == 1) {
        features <- match.arg(features, c("type", "state", "none"))
        features <- get(paste0(features, "_markers"))(x)
        if (length(features) == 0)
            stop("No features matched the specified marker class.")
    } else {
        # replace problematic characters
        features <- gsub("-", "_", features)
        features <- gsub(":", ".", features)
        stopifnot(features %in% rownames(x))
    }
    
    # calculate NRS
    cs_by_s <- split(seq_len(ncol(x)), x$sample_id)
    scores <- t(vapply(cs_by_s, function(cs)
        .nrs(assay(x, "exprs")[features, cs, drop = FALSE]),
        numeric(length(features))))
    
    # plot NRS in decreasing order
    md <- metadata(x)$experiment_info
    m <- match(levels(x$sample_id), md$sample_id)
    df <- melt(data.frame(scores, md[m, ]), id.var=names(md), 
        variable.name="antigen", value.name="NRS") 
    mean_scores <- colMeans(scores, na.rm=TRUE)
    o <- names(sort(mean_scores, decreasing=TRUE))
    df$antigen <- factor(df$antigen, levels=o)
    
    ggplot(df, aes_string(x="antigen", y="NRS")) +
        geom_point(aes_string(color=color_by), alpha=.8,
            position=position_jitter(width=.2, height=0)) +
        guides(color=guide_legend(override.aes=list(alpha=1))) +
        stat_summary(fun.y="mean", geom="point", 
            fill="white", size=2, shape=21) +
        geom_boxplot(width=.8, fill=NA, outlier.color=NA) +
        theme_classic() + theme(
            panel.grid.major=element_line(color="lightgrey", size=.25),
            axis.text=element_text(color="black"),
            axis.text.x=element_text(angle=90, vjust=.5, hjust=1))
}