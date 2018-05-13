# ==============================================================================
# Barplot of relative population abundances across samples & clusters
# ------------------------------------------------------------------------------
#' @rdname plotAbundances
#' @title Population frequencies across samples & clusters
#' 
#' @description 
#' Plots the relative population abundances of the specified clustering.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param k 
#'   specifies which clustering to use.
#' @param by 
#'   a character string specifying whether to plot 
#'   frequencies by samples or clusters.
#' @param group 
#'   a character string. Should corresponds to a column name of 
#'   \code{rowData(x)} other than "sample_id" and "cluster_id". 
#'   The default NULL will use the first factor available.
#' 
#' @return a \code{ggplot} object.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' # run clustering
#' lineage <- c("CD3", "CD45", "CD4", "CD20", "CD33", 
#'     "CD123", "CD14", "IgM", "HLA-DR", "CD7")
#' re <- cluster(re, cols_to_use=lineage)
#' # plot relative population abundances 
#' plotAbundances(re, k=12)                 # ...by sample 
#' plotAbundances(re, k=8, by="cluster_id") # ...by cluster
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importMethodsFrom S4Vectors metadata
# ------------------------------------------------------------------------------

setMethod(f="plotAbundances", 
    signature=signature(x="daFrame"), 
    definition=function(x, k=20, by=c("sample_id", "cluster_id"), group=NULL) {
    
        # validity checks
        by <- match.arg(by)
        k <- check_validity_of_k(x, k)
        valid <- setdiff(colnames(rowData(x)), c("sample_id", "cluster_id"))
        if (length(valid) == 0)
            stop("No factors to group by. Metadata should contain\n", 
                "at least one column other than 'file' and 'id'.")
        if (is.null(group)) {
            group <- valid[1]
        } else if (!is.character(group) | !group %in% valid) {
            stop("Argument 'group = ", dQuote(group), "' invalid.\n",
                "Should be one of: ", paste(dQuote(valid), collapse=", "))
        }
        
        # get cluster IDs & abundances
        cluster_ids <- cluster_codes(x)[, k][cluster_ids(x)]
        counts <- table(cluster_ids, sample_ids(x))
        
        # get frequencies by cluster & sample
        df <- melt(t(t(counts)/colSums(counts))*100, 
            varnames=c("cluster_id", "sample_id"),
            value.name="freq")
        
        # add metadata
        md <- metadata(x)$experiment_info
        m <- match(df$sample_id, md$sample_id)
        df <- data.frame(df, md[m, setdiff(names(md), names(df))])
        
        p <- ggplot(df, aes_string(y="freq")) +
            labs(x=NULL, y="Proportion [%]") + theme_bw() + theme(
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
                strip.background=element_rect(fill=NA, color=NA),
                strip.text=element_text(face="bold"),
                axis.ticks.x=element_blank(),
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
        
        switch(by,
            sample_id = p + facet_wrap(group, scales="free_x") +
                geom_bar(aes_string(x="sample_id", fill="factor(cluster_id)"), 
                    position="fill", stat="identity") +
                scale_fill_manual("cluster_id", values=cluster_cols) +
                scale_y_continuous(expand=c(0,0), labels=seq(0,100,25)) +
                theme(panel.border=element_blank()),
            cluster_id = p + facet_wrap(~cluster_id, scales="free_y", ncol=4) +
                guides(fill=FALSE) + geom_boxplot(aes_string(
                    x=group, color=group, fill=group),
                    position=position_dodge(), alpha=.25, outlier.color=NA) + 
                geom_point(position=position_jitter(width=.25),
                    aes_string(x=group, y="freq", color=group)) +
                theme(panel.grid.major=element_line(color="grey", size=.25))
        )
    }
)
