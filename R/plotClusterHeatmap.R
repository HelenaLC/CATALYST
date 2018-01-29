# ==============================================================================
# Heatmap of median marker expressions across clusters
# ------------------------------------------------------------------------------
#' @rdname plotClusterHeatmap
#' @title Median marker expressions across clusters
#' 
#' @description 
#' Plots a heatmap of median lineage marker expressions across clusters, 
#' and median functional marker expressions within each cluster across samples.
#'
#' @param x expression matrix.
#' @param scaled logical specifying whether scaled values should be plotted
#' (see below for details).
#' @param k specifies the clustering across which 
#' median marker expressions should be computed.
#' @param m specifies the second layer of clustering to be shown. 
#' @param type2 a character vector or \code{"all"}. 
#' Specifies for which type 2 markers heatmaps should be plotted.
#' @param out_path a character string. If specified, 
#' output will be generated in this location. Defaults to NULL.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @details Scaled values corresponds to cofactor arcsinh-transformed 
#' expression values scaled between 0 and 1 using 1% and 99% percentiles as 
#' boundaries. Hierarchical clustering is performed on the unscaled data.
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
#' plotClusterHeatmap(re, k=20, m=12, type2="pS6")
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import ComplexHeatmap
#' @importFrom dplyr funs group_by summarise_all
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 dcast
#' @importFrom stats dist
#' @export
# ==============================================================================

setMethod(f="plotClusterHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, scaled=TRUE, k=20, m=NULL, 
        type2="all", out_path=NULL) {
        
        check_validity_of_k(x, k)
        k <- as.character(k)
        if (!is.null(m)) {
            check_validity_of_k(x, m)
            m <- as.character(m)
        }
        es <- exprs(x)
        cluster_ids <- cluster_codes(x)[, k][cluster_ids(x)]
        n_clusters <- length(unique(cluster_codes(x)[, k]))
        
        # compute medians across clusters
        med_exprs <- data.frame(es, cluster_id=cluster_ids) %>%
            group_by(cluster_id) %>% summarize_all(funs(median))
        
        # cluster based on markers used for clustering
        d <- stats::dist(med_exprs[, type1(x)], method="euclidean")
        row_clustering <- hclust(d, method="average")
        
        # row labels and heatmap annotations
        row_anno <- factor(seq_len(n_clusters))
        row_cols <- cluster_cols[seq_len(nlevels(row_anno))]
        names(row_cols) <- levels(row_anno)
        heat_anno <- Heatmap(row_anno, row_cols, "cluster_id",
            cluster_rows=row_clustering, cluster_columns=FALSE,
            row_dend_reorder=FALSE, width=unit(.5, "cm"))
        
        # merging annotation
        merging_anno <- NULL
        if (!is.null(m)) {
            merging_ids <- factor(cluster_codes(x)[, m])[
                match(seq_len(n_clusters), cluster_codes(x)[, k])]
            merging_cols <- cluster_cols[seq_len(nlevels(merging_ids))]
            names(merging_cols) <- levels(merging_ids)
            merging_anno <- Heatmap(merging_ids, merging_cols, "merging_id",
                cluster_rows=row_clustering, cluster_columns=FALSE,
                row_dend_reorder=FALSE, width=unit(.5, "cm"))
        }
        
        # heatmap for type1 markers
        hm_cols <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100)
        if (scaled) {
            es0 <- scale_exprs(es)
            hm1_exprs <- data.frame(es0, cluster_id=cluster_ids) %>%
                group_by(cluster_id) %>% summarize_all(funs(median))
            hm2_exprs <- es0
        } else {
            hm1_exprs <- med_exprs
            hm2_exprs <- es
        }
        hm_l <- Heatmap(
            hm1_exprs[, type1(x)], hm_cols, "expression", 
            column_names_gp=gpar(fontsize=8),
            cluster_rows=row_clustering, cluster_columns=FALSE,
            heatmap_legend_param=list(color_bar="continuous"))
        
        # compute cluster frequencies
        counts <- as.numeric(table(cluster_ids))
        freqs <- round(counts/sum(counts)*100, 2)
        freq_labs <- paste0(seq_len(n_clusters), " (", freqs, "%)")
        freq_bars <- rowAnnotation(width=unit(2, "cm"), "Frequency [%]"=
                row_anno_barplot(x=freqs, border=FALSE, axis=TRUE,
                    gp=gpar(fill="grey50", col="white"), bar_with=.75))
        freq_anno <- rowAnnotation(
            text=row_anno_text(freq_labs), 
            width=max_text_width(freq_labs))
        
        p <- merging_anno + heat_anno + hm_l + freq_bars + freq_anno
        
        # heatmaps for type2 markers
        df <- data.frame(hm2_exprs[, type2(x)],
            sample_id=sample_ids(x), cluster_id=cluster_ids) %>%
            group_by(sample_id, cluster_id) %>% summarise_all(funs(median))
        if (type2 == "all") t2 <- type2(x) else t2 <- type2
        for (i in t2) {
            mat <- dcast(hm_exprs[, c("sample_id", "cluster_id", i)], 
                cluster_id~sample_id, value.var=i)[, -1]
            p <- p + Heatmap(mat, hm_cols, show_heatmap_legend=FALSE,
                column_title=i, column_names_gp=gpar(fontsize=8), 
                cluster_rows=row_clustering, cluster_columns=FALSE)
        }
        if (!is.null(out_path)) {
            out_nm <- file.path(out_path, "clustering_heatmap.pdf")
            pdf(out_nm, width=36, height=6); draw(p); dev.off()
        } else {
            draw(p)
        }
    }
)
