# ==============================================================================
# Barplot of event counts per sample
# ------------------------------------------------------------------------------

#' @rdname plotCounts
#' @title Number of events per sample
#' 
#' @description 
#' Barplot of the number of cells measured for each sample.
#'
#' @param x a \code{\link{daFrame}}.
#' @param color_by a character string that appears as a column name in the
#' metadata-table of the input \code{daFrame}. Specifies the color coding.
#' In the case of multiple conditions, \code{"condition"} will color points
#' according to each sample's combined conditions.
#' @param out_path a character string. If specified, 
#' output will be generated in this location. Defaults to NULL.
#' 
#' @return a \code{\link{ggplot}} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotCounts(re)
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import ggplot2 SummarizedExperiment
# ==============================================================================

setMethod(f="plotCounts", 
    signature=signature(x="daFrame"), 
    definition=function(x, color_by="condition") {
        md <- metadata(x)[[1]]
        n_events <- metadata(x)$n_events
        max <- ceiling(max(n_events)/1e3)*1e3 + 1e3
        conds <- grep("condition", colnames(md), value=TRUE)
        conds_combined <- apply(md[, conds], 1, paste, collapse="/")
        df <- data.frame(n_events, md, condition=conds_combined)
        ggplot(df, aes_string(x="sample_id", y="n_events",
            fill=color_by)) + geom_bar(stat="identity", width=.75) +  
            geom_label(aes_string(label="n_events", vjust=-.1), 
                fontface="bold.italic", fill=NA, label.size=0) +
            scale_y_continuous(limits=c(0,max), expand=c(0,0)) +
            theme_minimal() + theme(
                panel.grid.minor=element_blank(), 
                panel.grid.major.x=element_blank(),
                panel.grid.major.y=element_line(color="grey", size=.3), 
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, hjust=1, vjust=1))
    }
)
