# ==============================================================================
# Heatmap of median marker expressions across samples
# ------------------------------------------------------------------------------
#' @rdname plotMedExprs
#' @title Median marker expressions across samples
#' 
#' @description Plots median marker expressions across samples
#' computed on arcsinh-transformed intensities.
#'
#' @param x a \code{\link{daFrame}}.
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotMedExprs(re)
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @importFrom dplyr group_by summarize_all
#' @importFrom reshape2 dcast melt
# ==============================================================================

setMethod(f="plotMedExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x) {

        cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
        med_exprs <- data.frame(
            exprs(x)[, functional(x)],
            sample_id=sample_ids(x), 
            cluster_id=cluster_ids)
        med_exprs <- df %>% 
            group_by(sample_id, cluster_id) %>% 
            summarize_all(funs(median))
        med_exprs <- melt(med_exprs, 
            id.vars=c("sample_id", "cluster_id"),
            variable.name="antigen", 
            value.name="med_expr")
        med_exprs <- dcast(med_exprs, 
            cluster_id+antigen~sample_id, 
            value.var="med_expr")
        
        
        ggplot(df, aes_string(x="antigen", y="med_expr", col="condition")) +
            geom_point(size=2.5, alpha=.75, 
                position=position_jitter(width=.25, height=0)) +
            geom_boxplot(width=.75, fill=NA, outlier.color=NA) +
            stat_summary(fun.y="mean", geom="point", 
                fill="white", size=3, shape=21) +
            guides(color=guide_legend(override.aes=list(alpha=1))) +
            theme_classic() + theme(
                panel.grid.major=element_line(color="lightgrey", size=.25),
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, vjust=.5, hjust=1))
        
        
    }
)