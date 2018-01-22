# ==============================================================================
# Heatmap of median marker expressions across clusters
# ------------------------------------------------------------------------------
#' @rdname plotHeatmap
#' @title Median marker expressions across clusters
#' 
#' @description
#'
#' @param x expression matrix.
#' @param clustering numeric specifying the number of clusters across which 
#' median marker expression should be computed. If NULL, 
#' medians will be computed across samples. Defaults to NULL.
#' @param merging numeric specifying a second layer of clustering to be shown.
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
#' @import pheatmap
#' @importFrom dplyr funs group_by summarize_all
#' @importFrom RColorBrewer brewer.pal
#' @export
# ==============================================================================

setMethod(f="plotHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, clustering=NULL, merging=20, anno=FALSE) {
        
        # scale expressions for visualization
        es <- exprs(x)
        es0 <- scale_exprs(es)
    
        if (is.null(clustering)) {
            # compute medians across samples
            md <- metadata(x)[[1]]
            med_exprs <- data.frame(es, sample_id=sample_ids(x)) %>%
                group_by(sample_id) %>% summarize_all(funs(median))
            med_exprs <- t(data.frame(med_exprs, row.names=1))
            # get column annotations
            m <- match(colnames(med_exprs), md$sample_id)
            col_anno <- data.frame(row.names=colnames(med_exprs),
                condition=md$condition[m])
            
            pheatmap(med_exprs, annotation_col=col_anno,
                color=colorRampPalette(brewer.pal(8, "YlGnBu"))(100),
                display_numbers=anno, number_col="black")
        } else {
            # compute medians across clusters
            k_clustering <- paste0("k", clustering)
            cluster_ids <- cluster_ids(x)[, k_clustering]
            med_exprs <- data.frame(es, cluster_id=cluster_ids) %>%
                group_by(cluster_id) %>% summarize_all(funs(median))
            med_exprs_scaled <- data.frame(es0, cluster_id=cluster_ids) %>%
                group_by(cluster_id) %>% summarize_all(funs(median))
            
            # compute cluster frequencies
            counts <- table(cluster_ids)
            # cluster based on markers used for clustering
            d <- dist(med_exprs[, colnames(x)], method="euclidean")
            row_clustering <- hclust(d, method="average")
            # use scaled median expressions for plotting
            m <- as.matrix(med_exprs_scaled[, colnames(x)])
            rownames(m) <- med_exprs$cluster_id
            # row labels and heatmap annotations
            props <- paste0("(", round(counts/sum(counts)*100, 2), "%)")
            row_labs <- paste(rownames(m), props)
            row_anno <- merging_cols <- NULL
            
            if (!is.null(merging)) {
                inds <- match(seq_len(clustering), cluster_codes(x)[, k_clustering])
                new_ids <- cluster_codes(x)[, paste0("k", merging)][inds]
                new_ids <- factor(new_ids)
                row_anno <- data.frame(Merging=new_ids)
                n <- nlevels(new_ids)
                merging_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(n)
                names(merging_cols) <- levels(new_ids)
                cols_anno <- list(Merging=merging_cols)
            }
            cols <- c("darkslateblue", rev(brewer.pal(9, "RdYlBu")), "red3")
            pheatmap(m, color=colorRampPalette(cols)(100),
                border_color="white", cluster_colors=merging_cols,
                annotation_row=row_anno, annotation_color=cols_anno,
                cluster_rows=row_clustering, cluster_cols=FALSE,
                labels_row=row_labs, display_numbers=anno,
                number_color="black")
        }
    }
)
