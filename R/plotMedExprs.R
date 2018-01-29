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
#' 
#' # run clustering
#' lineage <- c("CD3", "CD45", "CD4", "CD20", "CD33", 
#'     "CD123", "CD14", "IgM", "HLA_DR", "CD7")
#' re <- cluster(re, cols_to_use=lineage)
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
    definition=function(x, k=20, type2="all", color_by="condition") {

        check_validity_of_k(x, k)
        
        if (type2 == "all") t2 <- type2(x) else t2 <- type2
        cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
        df <- data.frame(exprs(x)[, t2],
            sample_id=sample_ids(x), 
            cluster_id=cluster_ids)
        colnames(df)[seq_len(length(t2))] <- t2

        # compute medians across samples & clusters
        med_exprs <- df %>% 
            group_by(sample_id, cluster_id) %>% 
            summarize_all(funs(median))
        med_exprs <- melt(med_exprs, 
            id.vars=c("sample_id", "cluster_id"),
            variable.name="antigen", 
            value.name="med_expr")
        
        # get conditions
        md <- metadata(x)[[1]]
        m <- match(med_exprs$sample_id, md$sample_id)
        conds <- grep("condition", colnames(md), value=TRUE)
        for (i in seq_along(conds))
            med_exprs[, conds[i]] <- md[m, conds[i]]
        if (length(conds) > 1) {
            conds_combined <- apply(md[, conds], 1, paste, collapse="/")
            med_exprs$condition <- conds_combined[m]
        } 

        ggplot(med_exprs, aes_string(x="antigen", y="med_expr", col=color_by)) +
            facet_wrap(facets="cluster_id", scale="free_y") +
            geom_point(alpha=.5, 
                position=position_jitter(width=.25, height=0)) +
            geom_boxplot(width=.75, fill=NA, outlier.color=NA) +
            stat_summary(fun.y="mean", geom="point", 
                fill="white", size=2.5, shape=21) +
            guides(color=guide_legend(override.aes=list(alpha=1))) +
            ylab("median expression") + theme_classic() + theme(
                panel.grid.major=element_line(color="lightgrey", size=.25),
                strip.background=element_rect(fill=NA, color=NA),
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, vjust=.5, hjust=1))        
        
        ggplot(med_exprs, aes_string(x="cluster_id", y="med_expr", col=color_by)) +
            facet_wrap(facets="antigen", scale="free_y") +
            geom_point(alpha=.5, 
                position=position_jitter(width=.25, height=0)) +
            geom_boxplot(width=.75, fill=NA, outlier.color=NA) +
            stat_summary(fun.y="mean", geom="point", 
                fill="white", size=2.5, shape=21) +
            guides(color=guide_legend(override.aes=list(alpha=1))) +
            ylab("Median expression") + theme_classic() + theme(
                panel.grid.major=element_line(color="lightgrey", size=.25),
                strip.background=element_rect(fill=NA, color=NA),
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, vjust=.5, hjust=1))
    }
)