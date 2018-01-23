# ==============================================================================
# Barplot of relative population abundances across samples & clusters
# ------------------------------------------------------------------------------

#' @rdname plotAbundances
#' @title Population frequencies across samples & clusters
#' 
#' @description
#'
#' @param x a \code{\link{daFrame}}.
#' @param k specifies which clustering to use.
#' @param by a character string specifying whether 
#' to plot frequencies across samples or clusters.
#' 
#' @return a \code{ggplot} object.
#'
#' @references 
#' Nowicka M, Krieg C, Weber LM et al.
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md, merging_table)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' # plot relative population abundances 
#' plotAbundances(re, k=12)                 # ...across samples 
#' plotAbundances(re, k=8, by="cluster_id") # ...across clusters
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 SummarizedExperiment
# ==============================================================================

setMethod(f="plotAbundances", 
    signature=signature(x="daFrame"), 
    definition=function(x, k=20, by=c("sample_id", "cluster_id")) {
    
        cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
        counts <- table(cluster_ids, sample_ids(x))
        df <- data.frame(t(t(counts)/colSums(counts))*100)
        colnames(df) <- c("cluster_id", "sample_id", "freq")
        md <- metadata(x)[[1]]
        df$condition <- factor(md$condition[match(df$sample_id, md$sample_id)])
        
        p <- ggplot(df, aes_string(y="freq")) +
            labs(x=NULL, y="Proportion [%]") + theme_bw() + theme(
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
                strip.background=element_rect(fill=NA, color=NA),
                axis.ticks.x=element_blank(),
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
        switch(match.arg(by),
            sample_id = p + facet_wrap(~condition, scales="free_x") +
                geom_bar(aes_string(x="sample_id", fill="cluster_id"), 
                    position="fill", stat="identity") +
                scale_fill_manual(values=cluster_cols) +
                scale_y_continuous(expand=c(0,0), labels=seq(0,100,25)) +
                theme(panel.border=element_blank()),
            cluster_id = p + facet_wrap(~cluster_id, scales="free", ncol=4) +
                guides(fill=FALSE) + geom_boxplot(aes_string(
                    x="condition", color="condition", fill="condition"),
                    position=position_dodge(), alpha=.5, outlier.color=NA) + 
                geom_point(alpha=.75, position=position_jitter(),
                    aes_string(x="condition", y="freq", color="condition")) +
                theme(panel.grid.major=element_line(color="grey", size=.25))
        )
    }
)
