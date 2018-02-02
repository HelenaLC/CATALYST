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
#' @param hm2 a character string that specifies the right-hand side heatmap. 
#' Depending on the argument, the output will be as follows:
#' \itemize{
#' \item \code{"abundances"}: cluster frequencies across samples
#' \item \code{"type2"}: median type 2 marker expressions across clusters 
#' (analogous to the left-hand side heatmap)
#' \item a character string/vector corresponding to one/multiple marker(s): 
#' median marker expressions across samples and clusters}
#' @param k specifies the clustering across which 
#' median marker expressions should be computed.
#' @param m specifies the second layer of clustering to be shown. 
#' @param freq_bars,freq_labs logical. Should marginal histograms
#' of cluster frequencies be shown and annotated?
#' @param scaled logical specifying whether scaled values should be plotted
#' (see below for details).
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
#' plotClusterHeatmap(re, hm2="abundances")
#' plotClusterHeatmap(re, hm2="type2", k=12)
#' plotClusterHeatmap(re, hm2="pS6", k=12, m=8)
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
    definition=function(x, hm2=NULL, k=20, m=NULL, 
        freq_bars=TRUE, freq_labs=TRUE,
        scaled=TRUE, out_path=NULL) {
        
        check_validity_of_k(x, k)
        if (!is.null(hm2) && !hm2 %in% c("abundances", "type2", colnames(x)))
            stop("Invalid argument for 'hm2'. Should be NULL, ", 
                paste(dQuote(c("abundances", "type2")), collapse=", "), 
                " or a character string in ", 
                paste0("'colnames(", deparse(substitute(x)), ")'."))
        
        k <- as.character(k)
        if (!is.null(m)) {
            check_validity_of_k(x, m)
            m <- as.character(m)
        }
        
        es <- exprs(x)
        cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
        n_clusters <- nlevels(cluster_ids)
        
        # compute medians across clusters
        med_exprs <- data.frame(es, cluster_id=cluster_ids) %>%
            group_by(cluster_id) %>% summarize_all(funs(median))
        
        # cluster based on markers used for clustering
        d <- stats::dist(med_exprs[, type1(x)], method="euclidean")
        row_clustering <- hclust(d, method="average")
        
        # row labels and cluster annotations
        row_anno <- levels(cluster_ids)
        row_cols <- cluster_cols[seq_len(n_clusters)]
        names(row_cols) <- row_anno
        cluster_anno <- Heatmap(row_anno, row_cols, "cluster_id",
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
        
        # left-hand side heatmap:
        # type 1 median marker expressions across clusters
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
        hm1 <- Heatmap(hm1_exprs[, type1(x)], hm_cols, 
            paste0("scaled\n"[TRUE], "expression"), 
            column_names_gp=gpar(fontsize=8),
            cluster_rows=row_clustering, cluster_columns=FALSE,
            heatmap_legend_param=list(color_bar="continuous"))
        
        # compute cluster frequencies
        counts <- as.numeric(table(cluster_ids))
        freqs <- round(counts/sum(counts)*100, 2)
        if (freq_bars) {
            freq_bars <- rowAnnotation(width=unit(2, "cm"), "Frequency [%]"=
                    row_anno_barplot(x=freqs, border=FALSE, axis=TRUE,
                        gp=gpar(fill="grey50", col="white"), bar_with=.75))
        } else {
            freq_bars <- NULL
        }
        if (freq_labs) {
            freq_labs <- paste0(row_anno, " (", freqs, "%)")
            freq_anno <- rowAnnotation(
                text=row_anno_text(freq_labs), 
                width=max_text_width(freq_labs))
        } else {
            freq_anno <- NULL
        }

        p <- merging_anno + cluster_anno + hm1 + freq_bars + freq_anno
        
        # right-hand side heatmap
        if (!is.null(hm2)) {
            if (hm2 == "abundances") {
                # cluster frequencies across samples
                counts <- as.data.frame.matrix(
                    table(cluster_ids, sample_ids(x)))
                n_events <- colSums(counts)
                freqs <- t(t(counts) / n_events)
                p <- p + Heatmap(freqs, 
                    rev(brewer.pal(11, "PuOr")[2:10]), "frequency",
                    show_row_names=FALSE, column_names_gp=gpar(fontsize=8), 
                    cluster_rows=row_clustering, cluster_columns=FALSE,
                    heatmap_legend_param=list(color_bar="continuous"))
            } else if (hm2 == "type2") {
                # type 2 median expressions across clusters
                p <- p + Heatmap(hm1_exprs[, type2(x)], hm_cols, "", 
                    cluster_rows=row_clustering, cluster_columns=FALSE,
                    column_names_gp=gpar(fontsize=8), 
                    show_heatmap_legend=FALSE)
            } else {
                # median marker expressions across samples & clusters
                df <- data.frame(hm2_exprs,
                    sample_id=sample_ids(x), cluster_id=cluster_ids) %>%
                    group_by(sample_id, cluster_id) %>% 
                    summarise_all(funs(median))
                for (i in hm2) {
                    mat <- dcast(df[, c("sample_id", "cluster_id", i)], 
                        cluster_id~sample_id, value.var=i)[, -1]
                    p <- p + Heatmap(mat, hm_cols, show_heatmap_legend=FALSE,
                        column_title=i, column_names_gp=gpar(fontsize=8), 
                        cluster_rows=row_clustering, cluster_columns=FALSE)
                }
            }
        }
        
        if (!is.null(out_path)) {
            out_nm <- file.path(out_path, "clustering_heatmap.pdf")
            pdf(out_nm, width=36, height=6); draw(p); dev.off()
        } else {
            draw(p)
        }
    }
)
