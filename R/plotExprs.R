#' @rdname plotExprs
#' @title Plot expressions
#' 
#' @description 
#' Plots the smoothed densities of arcsinh-transformed marker intensities.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param color_by 
#'   character string. Has to appear as a column name of \code{rowData(x)}.
#'   Specifies the color coding.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return a \code{\link{ggplot}} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotExprs(re)
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
# ------------------------------------------------------------------------------

setMethod(f="plotExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x, color_by="condition") {
        
        df <- data.frame(exprs(x), rowData(x))
        df <- melt(df, id.vars=names(rowData(x)),
            variable.name="antigen", value.name="expression")
        
        ggplot(df, aes_string(x="expression", 
            col=color_by, group="sample_id"), fill=NULL) + 
            facet_wrap(~antigen, ncol=5, scales="free") +
            geom_density() + theme_classic() + theme(
                panel.grid=element_blank(), 
                strip.background=element_blank(),
                strip.text=element_text(face="bold"),
                axis.text=element_text(color="black"), 
                axis.title=element_text(color="black"))
    }
)

plot_clustering_distr_wrapper <- function(x, markers, k = 20, 
    clustering_distance="euclidean", clustering_linkage="average") {
    
    # validity checks
    if (!all(markers %in% colnames(x))) 
        stop("Invalid 'markers' specified.\n'",
            "all(markers %in% colnames(x)' should return TRUE.")
    k <- check_validity_of_k(x, k)
    
    # calculate median expression by cluster
    cluster_ids <- cluster_codes(x)[, k][cluster_ids(x)]
    df <- data.frame(exprs(x)[, markers], cluster_id=cluster_ids)
    med_exprs <- df %>% group_by_(~cluster_id) %>% summarize_all(funs(median))
    
    # sort clusters with hierarchical clustering
    d <- dist(med_exprs[, markers], method=clustering_distance)
    h <- hclust(d, method=clustering_linkage)
    o <- levels(cluster_ids)[rev(h$order)]
    
    # calculate cluster frequencies
    counts <- table(cluster_ids)
    freqs <- round(as.numeric(counts)/sum(counts)*100, 2)
    df$cluster_id <- factor(cluster_ids, 
        labels=paste0(levels(cluster_ids), " (", freqs, "%)")[as.numeric(o)])

    df_dat <- melt(df, id.vars="cluster_id", 
        value.name="expression", variable.name="antigen")
    df_dat$reference <- "no"
    # reference data
    df_ref <- df_dat
    df_ref$cluster_id <- "reference"
    df_ref$reference <- "yes"
    
    gg_df <- rbind(df_dat, df_ref)
    ggplot(gg_df, aes_string(x="expression", y="cluster_id", 
        col="reference", fill="reference")) + geom_density_ridges(alpha=.3) +
        facet_wrap(~antigen, scales="free_x", nrow=2) + theme_ridges() + 
        theme(legend.position="none", axis.title=element_text(size=10),
            axis.text=element_text(size=8), strip.text=element_text(size=8))
}