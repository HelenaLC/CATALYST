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
#' @return Plots the number of events measured per sample.
#'
#' @references 
#' Nowicka M, Krieg C, Weber LM et al.
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2
#' @importFrom flowCore fsApply
# ==============================================================================

setMethod(f="plotCounts", 
    signature=signature(x="daFrame"), 
    definition=function(x) {
        n_events <- fsApply(data(x), nrow)
        max <- ceiling(max(n_events)/1e3)*1e3 + 1e3
        df <- data.frame(sample_id=metadata(x)$sample_id, 
            n_events, condition=metadata(x)$condition)
        ggplot(df, aes_string(x="sample_id", y="n_events", fill="condition")) +
            geom_bar(stat="identity", width=.75) +  
            geom_label(aes_string(label="n_events"),
                label.size=0, fill="white", alpha=.5) +
            scale_y_continuous(limits=c(0,max), expand=c(0,0)) +
            theme_minimal() + theme(
                plot.title=element_text(face="bold"),
                panel.grid.minor=element_blank(), 
                panel.grid.major.x=element_blank(),
                panel.grid.major.y=element_line(color="black", size=.25), 
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, hjust=1, vjust=1))
    }
)