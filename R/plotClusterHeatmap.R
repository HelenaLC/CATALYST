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
#' @param split_by a character string that correpsonds to a column name
#' in the metadata table. If specified, the data will be subset accordingly
#' and multiple heatmaps will be plotted.
#' @param scale logical specifying whether scaled values should be plotted
#' (see below for details).
#' @param dend logical. Specifies if the row dendrogram should be drawn.
#' @param cluster_anno logical. Specifies if clusters should be annotated.
#' @param freq_bars,freq_labs logical. Should marginal histograms
#' of cluster frequencies be shown and annotated?
#' @param palette vector of colors to use. 
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
#' plotClusterHeatmap(re, hm2="type2", k=16, split_by='condition')
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
    definition=function(x, hm2=NULL, k=20, m=NULL, split_by=NULL, scale=TRUE, 
        dend=TRUE, cluster_anno=TRUE, freq_bars=TRUE, freq_labs=TRUE, 
        palette=rev(brewer.pal(11, "RdYlBu"))) {
        
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
        
        # specify color palette
        hm_cols <- colorRampPalette(palette)(100)
        
        es <- exprs(x)
        cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
        n_clusters <- nlevels(cluster_ids)
        
        # compute medians across clusters
        med_exprs <- data.frame(es, cluster_id=cluster_ids) %>%
            group_by(cluster_id) %>% summarize_all(funs(median))
        
        # cluster based on markers used for clustering
        d <- stats::dist(med_exprs[, type1(x)], method="euclidean")
        row_clustering <- hclust(d, method="average")
        
        # cluster & merging annotation
        if (cluster_anno) {
            row_anno <- levels(cluster_ids)
            row_cols <- cluster_cols[seq_len(n_clusters)]
            names(row_cols) <- row_anno
            cluster_anno <- Heatmap(
                matrix=row_anno, col=row_cols, name="cluster_id", 
                rect_gp=gpar(col="white"), width=unit(.5, "cm"),
                cluster_rows=row_clustering, cluster_columns=FALSE,
                show_row_dend=dend, row_dend_reorder=FALSE)
        }
        merging_anno <- NULL
        if (!is.null(m)) {
            merging_ids <- factor(cluster_codes(x)[, m])[
                match(seq_len(n_clusters), cluster_codes(x)[, k])]
            merging_cols <- cluster_cols[seq_len(nlevels(merging_ids))]
            names(merging_cols) <- levels(merging_ids)
            merging_anno <- Heatmap(
                matrix=merging_ids, col=merging_cols, name="merging_id", 
                rect_gp=gpar(col='white'), width=unit(.5, "cm"),
                cluster_rows=row_clustering, cluster_columns=FALSE,
                show_row_dend=dend, row_dend_reorder=FALSE)
        }
        
        # optional data subsetting
        if (is.null(split_by)) {
            many <- FALSE
            grps <- list(seq_len(nrow(x)))
        } else {
            # validity check 
            if (!split_by %in% colnames(rowData(x)))
                stop("Invalid argument 'split_by'.\nShould be one of ", 
                    paste(dQuote(colnames(rowData(x))), collapse=", "), " or NULL.")
            many <- TRUE
            grps <- split(seq_len(nrow(x)), rowData(x)[[split_by]])
        }
        
        hms <- sapply(seq_along(grps), function(i) {
            inds <- grps[[i]]
            # left-hand side heatmap:
            # type 1 median marker expressions across clusters
            if (scale) {
                es0 <- scale_exprs(es[inds, ])
                hm1_exprs <- data.frame(es0, cluster_id=cluster_ids[inds]) %>%
                    group_by(cluster_id) %>% summarize_all(funs(median))
                missing <- levels(cluster_ids)[!levels(cluster_ids) %in% hm1_exprs$cluster_id]
                hm1_exprs <- rbind(hm1_exprs, matrix(c(factor(missing), rep(rep(NA, ncol(es0)), length(missing))), 
                    nrow=length(missing), ncol=ncol(es0)+1, dimnames=list(NULL, colnames(hm1_exprs))))
                hm2_exprs <- es0
            } else {
                hm1_exprs <- med_exprs[inds, ]
                hm2_exprs <- es[inds, ]
            }
            hm1 <- Heatmap(matrix=hm1_exprs[, type1(x)], col=hm_cols, 
                name="expression", column_names_gp=gpar(fontsize=8),
                rect_gp=gpar(col='white'), na_col="lightgrey", 
                cluster_rows=row_clustering, cluster_columns=FALSE,
                heatmap_legend_param=list(color_bar="continuous"),
                show_row_dend=dend, column_title=names(grps)[i][many])
            
            # compute cluster frequencies
            bars <- anno <- NULL
            counts <- as.numeric(table(cluster_ids[inds]))
            freqs <- round(counts/sum(counts)*100, 2)
            if (freq_bars)
                bars <- rowAnnotation(width=unit(2, "cm"), "Frequency [%]"=
                        row_anno_barplot(x=freqs, border=FALSE, axis=TRUE,
                            gp=gpar(fill="grey50", col="white"), bar_with=.75))
            if (freq_labs) {
                labs <- paste0(row_anno, " (", freqs, "%)")
                anno <- rowAnnotation(
                    text=row_anno_text(labs), 
                    width=max_text_width(labs))
            }
            
            if (is.null(m) && class(hm1) != 'Heatmap') {
                p <- hm1 + bars + anno
            } else {
                p <- merging_anno + cluster_anno + hm1 + bars + anno
            }
            
            # right-hand side heatmap
            if (!is.null(hm2)) {
                if (hm2 == "abundances") {
                    # cluster frequencies across samples
                    counts <- as.data.frame.matrix(
                        table(cluster_ids[inds], sample_ids(x)[inds]))
                    n_events <- colSums(counts)
                    freqs <- t(t(counts) / n_events)
                    remove <- apply(freqs, 2, function(i) any(is.na(i)))
                    freqs <- freqs[, !remove]
                    p <- p + Heatmap(freqs, na_col="lightgrey",
                        rev(brewer.pal(11, "PuOr")[2:10]), "frequency",
                        rect_gp=gpar(col='white'),
                        show_row_names=FALSE, column_names_gp=gpar(fontsize=8), 
                        cluster_rows=row_clustering, cluster_columns=FALSE,
                        heatmap_legend_param=list(color_bar="continuous"))
                } else if (hm2 == "type2") {
                    # type 2 median expressions across clusters
                    p <- p + Heatmap(hm1_exprs[, type2(x)], hm_cols, 
                        rect_gp=gpar(col='white'),
                        na_col="lightgrey", show_heatmap_legend=FALSE,
                        cluster_rows=row_clustering, cluster_columns=FALSE,
                        column_names_gp=gpar(fontsize=8))
                } else {
                    # median marker expressions across samples & clusters
                    df <- data.frame(hm2_exprs, 
                        sample_id=sample_ids(x)[inds], 
                        cluster_id=cluster_ids[inds]) %>%
                        group_by(sample_id, cluster_id) %>% 
                        summarise_all(funs(median))
                    for (ch in hm2) {
                        mat <- dcast(df[, c("sample_id", "cluster_id", ch)], 
                            cluster_id~sample_id, value.var=ch)[, -1]
                        p <- p + Heatmap(mat, hm_cols, na_col="lightgrey",
                            rect_gp=gpar(col='white'),
                            cluster_rows=row_clustering, cluster_columns=FALSE,
                            show_heatmap_legend=FALSE, column_title=ch, 
                            column_names_gp=gpar(fontsize=8))
                    }
                }
                return(p)
            }
        })
        for (i in seq_along(hms)) 
            draw(hms[[i]])
    }
)
