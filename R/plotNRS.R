# ==============================================================================
# Marker ranking based on the non-redundancy score
# ------------------------------------------------------------------------------
#' @rdname plotNRS
#' @title Non-redundancy scores across markers
#' 
#' @description Plots non-redundancy scores (NRS) in decreasing order.
#'
#' @param x a \code{\link{daFrame}}
#' @param color_by a character string that appears as a column name in the
#' metadata-table of the input \code{daFrame}. Specifies the color coding.
#' In the case of multiple conditions, \code{"condition"} will color points
#' according to each sample's combined conditions.
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotNRS(re)
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
# ==============================================================================

setMethod(f="plotNRS", 
    signature=signature(x="daFrame"), 
    definition=function(x, color_by="condition") {
        
        valid_cols <- colnames(rowData(x))
        if (!color_by %in% c("condition", valid_cols)) {
            txt <- dQuote(c(setdiff(valid_cols, "condition"), "condition"))
            stop("Argument 'color_by = ", dQuote(color_by), "' invalid.\n",
                "       Should be one of: ", paste(txt, collapse=", "))
        }
        md <- metadata(x)[[1]]
        sample_ids <- md$sample_id
        
        # calculate NRS
        scores <- t(sapply(sample_ids, function(i) 
            nrs(exprs(x)[sample_ids(x) == i, ])))
        rownames(scores) <- sample_ids
        mean_scores <- colMeans(scores, na.rm=TRUE)
        
        # plot NRS in decreasing order
        o <- names(sort(mean_scores, decreasing=TRUE))
        scores <- data.frame(scores)
        scores$sample_id <- sample_ids
        df <- melt(scores, id.var="sample_id", 
            variable.name="antigen", value.name="NRS") 
        df$antigen <- factor(df$antigen, levels=o)
        
        m <- match(df$sample_id, sample_ids)
        md <- metadata(x)[[1]]
        conds <- grep("condition", colnames(md), value=TRUE)
        conds_combined <- apply(md[, conds], 1, paste, collapse="/")
        for (i in seq_along(conds))
            df[, conds[i]] <- md[, conds[i]][m]
        df$condition <- conds_combined[m]
        
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
