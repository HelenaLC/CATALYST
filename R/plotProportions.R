# ==============================================================================
# Barplot of relative population abundances across samples & clusters
# ------------------------------------------------------------------------------

#' @rdname plotProportions
#' @title Population frequencies across samples & clusters
#' 
#' @description
#'
#' @param fs a \code{\link{flowSet}} holding all samples.
#' @param md a \code{data.frame} containing the following columns: \itemize{
#' \item \code{file_name}: names of the fcs files contained in the input flowSet
#' \item \code{sample_id}: a unique identifier for each sample
#' \item \code{condition}: sample description (e.g. healthy or diseased)}
#' @param cluster_ids numeric vector of cluster ids for each event in \code{fs} 
#' @param by character specifying whether to plot relative abundances 
#' across samples or clusters
#' 
#' @return Barplot of relative population abundances.
#'
#' @references 
#' Nowicka M, Krieg C, Weber LM et al.
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' 
#' @import ggplot2
#' @importFrom flowCore flowSet fsApply
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @export
#' 
# ==============================================================================

setMethod(f="plotProportions", 
    signature=signature(fs="flowSet", md="data.frame", cluster_ids="numeric"), 
    definition=function(fs, md, cluster_ids, by=c("sample", "cluster")) {
    
        n_events <- fsApply(fs, nrow)
        m <- match(rownames(n_events), md$file_name)
        sample_ids <- rep(md$sample_id, n_events[m])
        counts <- table(cluster_ids, sample_ids)
        df <- data.frame(t(t(counts)/colSums(counts))*100)
        colnames(df) <- c("cluster_id", "sample_id", "freq")
        df$condition <- factor(md$condition[match(df$sample_id, md$sample_id)])
        
        p <- ggplot(df, aes_string(y="freq")) +
            labs(x=NULL, y="Proportion [%]") + theme_bw() + theme(
                panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="grey", size=.25),
                strip.background=element_rect(fill=NA, color=NA),
                axis.ticks.x=element_blank(),
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
        switch(match.arg(by),
            "sample" = p + facet_wrap(~condition, scales="free_x") +
                geom_bar(aes_string(x="sample_id", fill="cluster_id"), 
                    position="fill", stat="identity") +
                scale_fill_manual(values=cluster_cols) +
                scale_y_continuous(expand=c(0,0), labels=seq(0,100,25)) +
                theme(panel.border=element_blank()),
            "cluster" = p + facet_wrap(~cluster_id, scales="free", ncol=4) +
                geom_boxplot(aes_string(x="condition", 
                    color="condition", fill="condition"),
                    position=position_dodge(), alpha=.5, outlier.color=NA) + 
                geom_point(alpha=.75, position=position_jitter(),
                    aes_string(x="condition", y="freq", color="condition")) +
                guides(fill=FALSE)
        )
    }
)