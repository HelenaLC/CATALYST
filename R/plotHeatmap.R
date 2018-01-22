# ==============================================================================
# Heatmap of median marker expressions
# ------------------------------------------------------------------------------
#' @rdname plotHeatmap
#' @title Median marker expressions
#' 
#' @description
#'
#' @param x expression matrix.
#' @param clustering a numeric vector of length 2. 
#' Specifies the number of clusters across which median marker expressions 
#' should be computed and the second layer of clustering to be shown. 
#' If NULL (the default), medians will be computed across samples.
#' @param out_path a character string. If specified, 
#' output will be generated in this location. Defaults to NULL.
#' 
#' @return
#'
#' @references 
#' Nowicka M, Krieg C, Weber LM et al.
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return Heatmap of median marker expressions
#' 
#' @details Displayed values corresponds to cofactor 5 arcsinh-transformed 
#' expression values scaled between 0 and 1 using 1% and 99% percentiles as 
#' boundaries. Hierarchical clustering is performed on the unscaled data.
#' 
#' @examples
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @import ComplexHeatmap
#' @importFrom dplyr funs group_by summarize_all
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 dcast
#' @export
# ==============================================================================

setMethod(f="plotHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, clustering=NULL, out_path=NULL) {
        
        # scale expressions for visualization
        es <- exprs(x)
        es0 <- scale_exprs(es)
    
        if (is.null(clustering)) {
            # compute medians across samples
            md <- metadata(x)[[1]]
            med_exprs <- data.frame(es, sample_id=sample_ids(x)) %>%
                group_by(sample_id) %>% summarize_all(funs(median))
            med_exprs <- data.frame(med_exprs, row.names=1)
            # get column annotations
            m <- match(rownames(med_exprs), md$sample_id)
            cond_labs <- data.frame(
                row.names=rownames(med_exprs),
                condition=md$condition[m])
            cond_cols <- setNames(
                c("#F8766D", "#00BFC4"), 
                levels(md$condition))
            row_anno <- rowAnnotation(cond_labs, 
                "condition", list(condition=cond_cols),
                width=unit(.5, "cm"), gp=gpar(col="white"))
            # heatmap of medians across antigens and samples
            heat_cols <- colorRampPalette(brewer.pal(8, "YlGnBu"))(100)
            hm <- Heatmap(med_exprs, heat_cols, "expression",
                heatmap_legend_param=list(color_bar="continuous"),
                cell_fun=function(j,i,x,y,w,h,col){
                    grid.text(sprintf("%.2f", med_exprs[i,j]),x,y)})
            hm + row_anno
        } else {
            k1 <- clustering[1]
            k2 <- clustering[2]
            cluster_ids <- metadata(x)$cluster_codes[, k1][cluster_ids(x)]
            # compute medians across clusters
            med_exprs <- data.frame(es, cluster_id=cluster_ids) %>%
                group_by(cluster_id) %>% summarize_all(funs(median))
            med_exprs_scaled <- data.frame(es0, cluster_id=cluster_ids) %>%
                group_by(cluster_id) %>% summarize_all(funs(median))
            
            # cluster based on markers used for clustering
            d <- dist(med_exprs[, lineage(x)], method="euclidean")
            row_clustering <- hclust(d, method="average")
            
            # row labels and heatmap annotations
            row_anno <- factor(seq_len(k1))
            row_cols <- cluster_cols[seq_len(nlevels(row_anno))]
            names(row_cols) <- levels(row_anno)
            heat_anno <- Heatmap(row_anno, row_cols, "cluster_id",
                cluster_rows=row_clustering, cluster_columns=FALSE,
                row_dend_reorder=FALSE, width=unit(.5, "cm"))
            
            # merging annotation
            m <- match(seq_len(k1), metadata(x)$cluster_codes[, k1])
            merging_ids <- factor(metadata(x)$cluster_codes[, k2][m])
            merging_cols <- cluster_cols[seq_len(nlevels(merging_ids))]
            names(merging_cols) <- levels(merging_ids)
            merging_anno <- Heatmap(merging_ids, merging_cols, "merging_id",
                cluster_rows=row_clustering, cluster_columns=FALSE,
                row_dend_reorder=FALSE, width=unit(.5, "cm"))
            
            # heatmap for lineage markers
            heat_cols <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100)
            hm_l <- Heatmap(med_exprs_scaled[, lineage(x)], heat_cols, "expression", 
                column_title="Lineage markers", column_names_gp=gpar(fontsize=8),
                cluster_rows=row_clustering, cluster_columns=FALSE,
                heatmap_legend_param=list(color_bar="continuous"))
            
            # compute cluster frequencies
            counts <- as.numeric(table(cluster_ids))
            freqs <- round(counts/sum(counts)*100, 2)
            freq_labs <- paste0(seq_len(k1), " (", freqs, "%)")
            freq_bars <- rowAnnotation(width=unit(2, "cm"), "Frequency [%]"=
                    row_anno_barplot(x=freqs, border=FALSE, axis=TRUE,
                        gp=gpar(fill="grey50", col="white"), bar_with=.75))
            freq_anno <- rowAnnotation(
                    text=row_anno_text(freq_labs), 
                    width=max_text_width(freq_labs))

            p <- merging_anno + heat_anno + hm_l + freq_bars + freq_anno
            
            # heatmaps for functional markers
            df <- data.frame(es0[, functional(x)], 
                sample_id=sample_ids(x), cluster_id=cluster_ids) %>%
                group_by(sample_id, cluster_id) %>% summarise_all(funs(median))
            for (i in functional(x)) {
                mat <- dcast(df[, c("sample_id", "cluster_id", i)], 
                    cluster_id~sample_id, value.var=i)[, -1]
                p <- p + Heatmap(mat, heat_cols, show_heatmap_legend=FALSE,
                    column_title=i, column_names_gp=gpar(fontsize=8), 
                    cluster_rows=row_clustering, cluster_columns=FALSE)
            }
            if (!is.null(out_path)) {
                out_nm <- file.path(out_path, "clustering_heatmap.pdf")
                pdf(out_nm, width=36, height=6); draw(p); dev.off()
            } else {
                darw(p)
            }
        }
    }
)
