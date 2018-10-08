#' @rdname plotNRS
#' @title Plot non-redundancy scores
#' 
#' @description 
#' Plots non-redundancy scores (NRS) by markers in decreasing order.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param markers
#'   character string specifying which markers to include. Defaults to NULL 
#'   (= all markers). Alternatively, if the \code{colData(x)$marker_class} 
#'   column is specified, can be one of "type", "state", or "none".
#' @param color_by 
#'   character string. Has to appeara as a column name of \code{rowData(x)}.
#'   Specifies the color coding.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
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
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotNRS(re)
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
# ------------------------------------------------------------------------------

setMethod(f="plotNRS", 
    signature=signature(x="daFrame"), 
    definition=function(x, markers=NULL, color_by="condition") {
        
        # validity check
        md <- metadata(x)$experiment_info
        valid <- names(md)
        if (!color_by %in% valid)
            stop("Argument 'color_by = ", dQuote(color_by), "' invalid.\n",
                "Should be one of: ", paste(dQuote(valid), collapse=", "))

        # check validity of argument 'markers'
        if (is.null(markers)) {
            markers <- colnames(x)
        } else if (length(markers) == 1 &&
                markers %in% levels(colData(x)$marker_class)) {
            idx <- colData(x)$marker_class == markers
            if (!any(idx))
                stop(sprintf("No markers matched marker class '%s'.", markers))
            markers <- colnames(x)[idx]
        } else {
            # replace problematic characters
            markers <- gsub("-", "_", markers)
            markers <- gsub(":", ".", markers)
            stopifnot(all(markers %in% colnames(x)))
        }
        
        # calculate NRS
        scores <- t(sapply(md$sample_id, function(i) 
            nrs(exprs(x)[sample_ids(x) == i, markers])))
        mean_scores <- colMeans(scores, na.rm=TRUE)
        
        # plot NRS in decreasing order
        o <- names(sort(mean_scores, decreasing=TRUE))
        m <- match(rownames(scores), md$sample_id)
        scores <- data.frame(scores[, o], md[m, ])
        df <- melt(scores, id.var=names(md), 
            variable.name="antigen", value.name="NRS") 
        df$antigen <- factor(df$antigen, levels=o)

        ggplot(df, aes_string(x="antigen", y="NRS")) +
            geom_point(aes_string(color=color_by), alpha=.75,
                position=position_jitter(width=.25, height=0)) +
            guides(color=guide_legend(override.aes=list(alpha=1))) +
            stat_summary(fun.y="mean", geom="point", 
                fill="white", size=2, shape=21) +
            geom_boxplot(width=.75, fill=NA, outlier.color=NA) +
            theme_classic() + theme(
                panel.grid.major=element_line(color="lightgrey", size=.25),
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, vjust=.5, hjust=1))
    }
)
