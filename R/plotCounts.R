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
    definition=function(x) {
        md <- metadata(x)[[1]]
        n_events <- metadata(x)$n_events
        max <- ceiling(max(n_events)/1e3)*1e3 + 1e3
        df <- data.frame(n_events, 
            sample_id=md$sample_id, 
            condition=md$condition)
        ggplot(df, aes_string(x="sample_id", y="n_events", fill="condition")) +
            geom_bar(stat="identity", width=.75) +  
            geom_label(aes_string(label="n_events", vjust=-.25), 
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
