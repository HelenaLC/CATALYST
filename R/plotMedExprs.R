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
#' @param k a numeric or charactering specifying the clustering to use.
#' If \code{facette = "antigen"} this argument will be ignored.
#' @param facette character string \code{"antigen"} or \code{"cluster_id"}. 
#' Note that the latter requires having run \code{\link{cluster}} first.
#' @param group_by character string specifying the condition samples should be 
#' grouped by. If there are multiple conditions specified in the metadata-table,
#' \code{"condition"} will group samples according to their combined conditions.
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # plot median expressions
#' plotMedExprs(re)
#' 
#' # run clustering
#' lineage <- c("CD3", "CD45", "CD4", "CD20", "CD33", 
#'     "CD123", "CD14", "IgM", "HLA_DR", "CD7")
#' re <- cluster(re, cols_to_use=lineage)
#' 
#' # plot median expressions across clusters
#' plotMedExprs(re, facette="cluster_id", k=8)
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
#' @importFrom reshape2 melt
# ==============================================================================

setMethod(f="plotMedExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x, k=20, 
        facette=c("antigen", "cluster_id"), group_by="condition") {

        facette <- match.arg(facette)
        
        style <- list(ylab("median expression"),
            guides(color=guide_legend(override.aes=list(alpha=1))),
            geom_point(alpha=.75, 
                position=position_jitter(width=.25, height=0)),
            geom_boxplot(alpha=.5, width=.75, fill=NA, outlier.color=NA),
            theme_bw(), theme(
                axis.text=element_text(color="black"),
                panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="lightgrey", size=.25),
                strip.background=element_rect(fill="grey90", color=NA)))
        
        if (facette == "antigen") {
            med_exprs <- data.frame(exprs(x), sample_id=sample_ids(x)) %>% 
                group_by(sample_id) %>% summarize_all(funs(median))
            med_exprs <- melt(med_exprs, id.vars=c("sample_id"),
                variable.name="antigen", value.name="med_expr")
        } else{
            check_validity_of_k(x, k)
            cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
            # compute medians across samples & clusters
            med_exprs <- data.frame(exprs(x)[, state_markers(x)], 
                sample_id=sample_ids(x), cluster_id=cluster_ids) %>% 
                group_by(sample_id, cluster_id) %>% summarize_all(funs(median))
            med_exprs <- melt(med_exprs, id.vars=c("sample_id", "cluster_id"),
                variable.name="antigen", value.name="med_expr")
        }
        # get conditions
        md <- metadata(x)[[1]]
        m <- match(med_exprs$sample_id, md$sample_id)
        conds <- grep("condition", colnames(md), value=TRUE)
        conds_df <- cbind(sapply(conds, function(i) factor(md[m, i])))
        med_exprs <- data.frame(med_exprs, conds_df)
        if (length(conds) > 1) {
            conds_combined <- apply(md[, conds], 1, paste, collapse="/")
            med_exprs$condition <- conds_combined[m]
        }
        
        if (facette == "antigen") {
            ggplot(med_exprs, 
                aes_string(x=group_by, y="med_expr", col=group_by)) +
                stat_summary(fun.y="mean", geom="point", size=2.5, shape=21) +
                facet_wrap(facets="antigen", scales="free_y") + style + theme(
                axis.text.x=element_blank(), axis.ticks.x=element_blank())
        } else {
            ggplot(med_exprs, 
                aes_string(x="antigen", y="med_expr", col=group_by)) +
                facet_wrap(facets="cluster_id", scales="free_y") + style + 
                theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
        }
    }
)
        