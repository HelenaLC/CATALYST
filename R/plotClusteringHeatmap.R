# ==============================================================================
# Heatmap of median marker expressions across clusters
# ------------------------------------------------------------------------------
#' @rdname plotClusteringHeatmap
#' @title Median marker expressions across clusters
#' 
#' @description
#'
#' @param fs 
#' @param cluster_ids
#' @param merging_ids
#' @param anno logical indicating whether to display values inside each bin
#' 
#' @return
#'
#' @references 
#' Nowicka M, Krieg C, Weber LM et al.
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return Heatmap of median marker expressions across clusters.
#' 
#' @details Displayed values corresponds to cofactor 5 arcsinh-transformed 
#' expression values scaled between 0 and 1 using 1% and 99% percentiles as 
#' boundaries. Hierarchical clustering is performed on the unscaled data.
#' 
#' @examples
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @import ggplot2
#' @importFrom flowCore flowSet fsApply
#' @importFrom matrixStats colQuantiles
#' @export
# ==============================================================================

setMethod(f="plotClusteringHeatmap", 
    signature=signature(x="matrix", cluster_ids="vector"), 
    definition=function(x, cluster_ids, merging_ids=NULL, anno=FALSE) {
        
        # validity check
        n_events <- nrow(x)
        n_ids <- length(cluster_ids)
        if (n_ids != n_events) 
            stop(cluster_ids, " cluster ID(s) provided but ", 
                n_events, " needed.", call.=FALSE)
        
        # scale expressions for visualization
        x0 <- scale_exprs(x)
    
        # compute medians across clusters
        med_exprs <- data.frame(x, cluster_id=cluster_ids) %>%
            group_by(cluster_id) %>% summarize_all(funs(median))
        med_exprs_scaled <- data.frame(x0, cluster_id=cluster_ids) %>%
            group_by(cluster_id) %>% summarize_all(funs(median))
        
        # compute cluster frequencies
        counts <- table(cluster_ids)
        # cluster based on markers used for clustering
        d <- dist(med_exprs[, colnames(x)], method="euclidean")
        row_clustering <- hclust(d, method="average")
        # use scaled median expressions for plotting
        m <- as.matrix(med_exprs_scaled[, colnames(x0)])
        rownames(m) <- med_exprs$cluster_id
        # row labels and heatmap annotations
        props <- paste0("(", round(counts/sum(counts)*100, 2), "%)")
        row_labs <- paste(rownames(m), props)
        row_anno <- merging_cols <- NULL
        if (!is.null(merging_ids)) {
            merging_ids <- factor(merging_ids)
            row_anno <- data.frame(Merging=merging_ids)
            n <- nlevels(merging_ids)
            merging_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(n)
            names(merging_cols) <- levels(merging_ids)
            cols_anno <- list(Merging=merging_cols)
        }
        cols <- c("darkslateblue", rev(brewer.pal(9, "RdYlBu")), "red3")
        pheatmap(m, main="Scaled median marker intensities", 
            color=colorRampPalette(cols)(100),
            border_color="white", cluster_colors=merging_cols,
            annotation_row=row_anno, annotation_color=cols_anno,
            cluster_rows=row_clustering, cluster_cols=FALSE,
            labels_row=row_labs, display_numbers=anno,
            number_color="black")
    }
)