#' @rdname plotClusterHeatmap
#' @title Plot cluster heatmap
#' 
#' @description 
#' Plots heatmaps summarizing a clustering and/or metaclustering of interest.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param hm2 
#'   character string. Specifies the right-hand side heatmap. One of:
#'   \itemize{
#'   \item{\code{"abundances"}: cluster frequencies across samples}
#'   \item{\code{"state_markers"}: median cell state marker expressions 
#'     across clusters (analogous to the left-hand side heatmap)}
#'   \item{a character string/vector corresponding to one/multiple marker(s): 
#'     median marker expressions across samples and clusters}}
#' @param k 
#'   character string. Specifies the clustering 
#'   across which median marker expressions should be computed.
#' @param m 
#'   character string. Specifies the metaclustering to be shown. 
#'   (This is for display only and will not effect any computations!) 
#' @param fun
#'   character string specifying the function to use as summary statistic.
#' @param cluster_anno 
#'   logical. Specifies if clusters should be annotated.
#' @param split_by 
#'   character string. Must corresponds to a column name of \code{rowData(x)}. 
#'   If specified, the data will be subset according to this variable, 
#'   and multiple heatmaps will be drawn.
#' @param scale 
#'   logical. Specifies whether scaled values should be plotted.
#'   (see below for details)
#' @param draw_dend 
#'   logical. Specifies if the row dendrogram should be drawn.
#' @param draw_freqs 
#'   logical. Specifyies whether to display cell counts and proportions.
#' @param palette 
#'   character vector of colors to interpolate. 
#' 
#' @details Scaled values corresponds to cofactor arcsinh-transformed 
#' expression values scaled between 0 and 1 using 1% and 99% percentiles as 
#' boundaries. Hierarchical clustering is performed on the unscaled data.
#' 
#' In its 1st panel, \code{plotClusterHeatmap} will display
#' median (scaled, arcsinh-transformed) cell-type marker expressions (across all samples).
#' Depending on argument \code{hm2}, the 2nd panel will contain one of:
#' \itemize{
#' \item{relataive cluster abundances by sample}
#' \item{median (scaled, arcsinh-transformed) 
#'   cell-state marker expressions (across all samples)}
#' \item{median (scaled, arcsinh-transformed) 
#'   cell-state marker expressions by sample}
#' }
#' 
#' @return a \code{\link{HeatmapList-class}} object.
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
#' 
#' # run clustering
#' re <- cluster(re)
#' 
#' plotClusterHeatmap(re, hm2="abundances")
#' plotClusterHeatmap(re, hm2="abundances", draw_freqs=TRUE)
#' plotClusterHeatmap(re, hm2="state_markers", k="meta16", split_by='condition')
#' plotClusterHeatmap(re, hm2="pS6", k="meta12", m="meta8")
#' plotClusterHeatmap(re, hm2="abundances", scale=FALSE, draw_freqs=TRUE)
#' 
#' @import ComplexHeatmap
#' @importFrom dplyr group_by_ summarise_all
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 acast
#' @importFrom S4Vectors metadata
#' @importFrom stats dist
# ------------------------------------------------------------------------------

setMethod(f="plotClusterHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, hm2=NULL, k="meta20", m=NULL, fun=c("median", "mean"),
        cluster_anno=TRUE, split_by=NULL, scale=TRUE, draw_dend=TRUE, draw_freqs=FALSE, 
        palette=rev(brewer.pal(11, "RdYlBu"))) {
        
        fun <- match.arg(fun)
        
        # check validity of argument 'hm2'
        if (!is.null(hm2) && !hm2 %in% 
                c("abundances", "state_markers", colnames(x)))
            stop("Invalid argument for 'hm2'. Should be NULL, ", 
                paste(dQuote(c("abundances", "state_markers")), 
                    collapse=", "), " or a character string in ", 
                paste0("'colnames(", deparse(substitute(x)), ")'."))
        
        # check validity of arguments 'k' and 'm'
        k <- check_validity_of_k(x, k)
        m <- check_validity_of_k(x, m)
        
        cluster_ids <- cluster_codes(x)[cluster_ids(x), k]
        n_clusters <- nlevels(cluster_ids)
        
        # medians marker exprs. across clusters
        med_exprs <- data.frame(exprs(x), cluster_ids) %>%
            group_by_(~cluster_ids) %>% summarize_all(fun)
        
        # hierarchical clustering on cell-type markers
        d <- stats::dist(med_exprs[, type_markers(x)])
        row_clustering <- stats::hclust(d, method="average")
        
        # clustering row annotation 
        if (cluster_anno) {
            anno <- levels(cluster_ids)
            if (n_clusters > 30) {
                cols <- colorRampPalette(cluster_cols)(n_clusters)
            } else {
                cols <- cluster_cols[seq_len(n_clusters)]
            }
            cols <- setNames(cols, anno)
            cluster_anno <- row_anno(anno, cols, 
                "cluster_id", row_clustering, draw_dend)
        }
        # merging row annotation
        if (length(m) != 0) {
            anno <- factor(cluster_codes(x)[, m][match(
                seq_len(n_clusters), cluster_codes(x)[, k])])
            if (nlevels(anno) > 30) {
                cols <- colorRampPalette(cluster_cols)(nlevels(anno))
            } else {
                cols <- cluster_cols[seq_len(nlevels(anno))]
            }
            cols <- setNames(cols, levels(anno))
            merging_anno <- row_anno(anno, cols, 
                "merging_id", row_clustering, draw_dend)
        }
        
        # subsetting
        if (is.null(split_by)) {
            many <- FALSE
            groups <- list(seq_len(nrow(x)))
        } else {
            # validity check
            valid_opts <- colnames(rowData(x))
            if (!split_by %in% valid_opts)
                stop("Invalid argument 'split_by'.\nShould be one of ", 
                    paste(dQuote(valid_opts), collapse=", "), " or NULL.")
            many <- TRUE
            groups <- split(seq_len(nrow(x)), rowData(x)[[split_by]])
        }
        
        hm_cols <- colorRampPalette(palette)(100)
        hms <- sapply(seq_along(groups), function(i) {
            inds <- groups[[i]]
            # left-hand side heatmap:
            # median cell-type marker expressions across clusters
            if (scale) {
                es0 <- scale_exprs(exprs(x)[inds, , drop=FALSE])
                hm1_es <- data.frame(es0, cluster_ids=cluster_ids[inds]) %>% 
                    group_by_(~cluster_ids) %>% summarize_all(fun)
                hm2_es <- es0
            } else {
                if (!many) {
                    hm1_es <- med_exprs
                    hm2_es <- exprs(x)
                } else {
                    hm2_es <- exprs(x)[inds, , drop=FALSE]
                    hm1_es <- data.frame(hm2_es, 
                        cluster_ids=cluster_ids[inds]) %>% 
                        group_by_(~cluster_ids) %>% summarize_all(fun)
                }
            }
            
            # add clusters if any missing
            missing <- levels(cluster_ids)[
                !levels(cluster_ids) %in% hm1_es$cluster_ids]
            if (length(missing) > 0) {
                na_matrix <- matrix(NA, 
                    nrow=length(missing), ncol=ncol(hm1_es)-1,
                    dimnames=list(NULL, colnames(hm1_es)[-1]))
                na_df <- data.frame(cluster_ids=missing, na_matrix)
                hm1_es <- rbind(hm1_es, na_df) %>% arrange(cluster_ids)
            }
            
            hm1 <- Heatmap(matrix=hm1_es[, type_markers(x)], col=hm_cols, 
                name="expression", column_names_gp=gpar(fontsize=8),
                rect_gp=gpar(col='white'), na_col="lightgrey", 
                cluster_rows=row_clustering, cluster_columns=FALSE,
                show_row_dend=draw_dend, column_title=names(groups)[i][many])
            
            # cluster frequencies
            freq_bars <- freq_anno <- NULL
            if (draw_freqs) {
                counts <- as.numeric(table(cluster_ids[inds]))
                freqs <- round(counts/sum(counts)*100, 2)
                freq_bars <- rowAnnotation("Frequency [%]"=row_anno_barplot(
                    x=freqs, axis=TRUE, border=FALSE, bar_with=.8, 
                    gp=gpar(fill="grey50", col="white")), width=unit(2, "cm"))
                labs <- paste0(levels(cluster_ids), " (", freqs, "%)")
                freq_anno <- rowAnnotation(
                    text=row_anno_text(labs), 
                    width=max_text_width(labs))
            }
            
            # combine row annotations, heatmap, 
            # and frequency bars & labels
            p <- hm1 + freq_bars + freq_anno
            if (is(cluster_anno, "Heatmap")) 
                p <- cluster_anno + p
            if (exists("merging_anno")) 
                p <- merging_anno + p
            
            # right-hand side heatmap
            if (!is.null(hm2)) {
                if (hm2 == "abundances") {
                    # cluster frequencies across samples
                    counts <- as.data.frame.matrix(table(
                        cluster_ids[inds], sample_ids(x)[inds]))
                    freqs <- t(t(counts) / colSums(counts))
                    keep <- !apply(freqs, 2, function(x) all(is.na(x)))
                    freqs <- matrix(
                        freqs[, keep, drop=FALSE], nlevels(cluster_ids[inds]), 
                        dimnames=list(NULL, names(keep)[keep]))
                    p <- p + Heatmap(matrix=freqs,
                        name="frequency", na_col="lightgrey", 
                        rect_gp=gpar(col="white"), 
                        show_row_names=FALSE, 
                        column_names_gp=gpar(fontsize=8), 
                        cluster_rows=row_clustering, cluster_columns=FALSE)
                } else if (hm2 == "state_markers") {
                    # median cell state marker expressions across clusters
                    p <- p + Heatmap(show_heatmap_legend=FALSE, 
                        matrix=hm1_es[, state_markers(x)], col=hm_cols, 
                        na_col="lightgrey", rect_gp=gpar(col='white'), 
                        cluster_rows=row_clustering, cluster_columns=FALSE,
                        column_names_gp=gpar(fontsize=8))
                } else {
                    # median marker expression across samples & clusters
                    meds <- data.frame(hm2_es, 
                        sample_id=sample_ids(x)[inds], 
                        cluster_id=cluster_ids[inds]) %>%
                        group_by_(~sample_id, ~cluster_id) %>% 
                        summarise_all(fun)
                    for (ch in hm2) {
                        ch_meds <- acast(
                            meds[, c("sample_id", "cluster_id", ch)], 
                            formula=cluster_id~sample_id, value.var=ch)
                        p <- p + Heatmap(matrix=ch_meds, col=hm_cols, 
                            na_col="lightgrey", rect_gp=gpar(col='white'),
                            show_heatmap_legend=FALSE, show_row_names=FALSE,
                            cluster_rows=row_clustering, cluster_columns=FALSE,
                            column_title=ch, column_names_gp=gpar(fontsize=8))
                    }
                }
            }
            return(p)
        })
        for (i in seq_along(hms)) 
            draw(hms[[i]])
        invisible(hms)
    }
)