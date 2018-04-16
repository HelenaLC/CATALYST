# ==============================================================================
# Heatmap for differental abundance & state analysis
# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
#' @title Median marker expressions across clusters
#' 
#' @param type character string. 
#' Specifies whether to plot a heatmap for differential abundance (DA) 
#' or differential state (DS) test results.
#' @param top_n numeric. 
#' Number of top clusters (if \code{type = "DA"}) or cluster-marker 
#' combinations (if \code{type = "DS"}) to display. Defaults to 20.
#' @param all logical. 
#' Specifies whether to display all clusters or cluster-marker combinations. 
#' If \code{TRUE}, \code{top_n} will be ignored.
#' @param order logical.
#' Specifies whether results should be ordered by adjusted p-value.
#' 
#' @author
#' Lukas M Weber, Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @import ComplexHeatmap
#' @import dplyr
#' @importFrom circlize colorRamp2
#' @importFrom stats quantile
#' @export
# ------------------------------------------------------------------------------

setMethod(f="plotDiffHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, dt_res, type=c("DA", "DS"), 
        top_n=20, all=FALSE, order=TRUE) {
        
        type <- match.arg(type, c("DA", "DS"))
        
        # get no. of clusters & cluster IDs
        k <- nlevels(rowData(dt_res)$cluster_id)
        cluster_ids <- cluster_codes(x)[, k][cluster_ids(x)]
        df <- data.frame(exprs(x), 
            sample_id=sample_ids(x), 
            cluster_id=cluster_ids)
        
        # compute medians across samples & clusters
        meds_across_samples <- df %>% group_by(sample_id) %>% 
            summarize_at(colnames(x), funs(median))
        meds_across_clusters <- df %>% group_by(cluster_id) %>% 
            summarise_at(colnames(x), funs(median))

        # color scale: 
        # 1%, 50%, 99% percentiles across medians & markers
        qs <- quantile(unlist(meds_across_clusters),
            c(0.01, 0.5, 0.99), TRUE)
        hm_cols <- colorRamp2(qs, c("royalblue3", "white", "tomato2"))
        
        # get clusters or cluster-marker combinations to plot
        dt <- rowData(dt)
        if (order) {
            cols <- colnames(dt) %in% c("FDR", "adj.P.Val", "p_adj")
            dt <- dt[order(dt[, cols]), , drop=FALSE]
        }
        if (all | top_n > nrow(dt)) 
            top_n <- nrow(dt)
        top <- dt[seq_len(top_n), ]
        meds_across_clusters <- meds_across_clusters[
            match(top$cluster_id, rownames(meds_across_clusters)), ]
        
        # 1st heatmap:
        # median cell-type marker expressions across clusters
        hm1 <- Heatmap(col=hm_cols, name="expression", 
            matrix=meds_across_clusters[, type_markers(x)],
            column_title="type_markers", column_title_side="bottom", 
            cluster_columns=FALSE, row_names_side="left", 
            clustering_distance_rows="euclidean",
            clustering_method_rows="median")
        
        meds <- df %>% group_by(cluster_id, sample_id) %>% summarise_all(funs(median))
        meds <- setNames(lapply(colnames(x), function(marker)
            reshape2::acast(meds, cluster_id~sample_id, 
                value.var=marker, fill=0)), colnames(x))
        meds <- meds[top$cluster_id, , drop=FALSE]
       
        # 2nd heatmap:
        # type = "DA": cluster sizes across samples
        # type = "DS": median cell-state marker expressions across clusters
        if (type == "DA") {
            counts <- df %>% count(cluster_id, sample_id) %>% complete(sample_id)
            counts <- reshape2::acast(counts, cluster_id~sample_id, value.var="n", fill=0)
            counts <- counts[top$cluster_id, , drop=FALSE]
            
            hm2 <- Heatmap(matrix=counts, name="n_cells",
                col=colorRamp2(range(counts), c("black", "gold")),
                cluster_columns=FALSE, show_row_names=FALSE)
        } else if (type == "DS") { 
            hm2 <- Heatmap(col=colors, show_heatmap_legend=FALSE,
                matrix=meds_across_clusters[, state_markers(x)],
                column_title="state_markers", 
                column_title_side="bottom", 
                cluster_columns=FALSE, 
                row_names_side = "left", 
                clustering_distance_rows="euclidean",
                clustering_method_rows="median")
        }
        
        # if type = "DS", 3rd heatmap:
        # median cell-state marker expressions across samples
        
        # subset top markers & clusters
        data <- mapply(function(marker, ids) marker[ids, , drop = FALSE], 
            meds[top$marker], top$cluster_id, SIMPLIFY=FALSE)
        data <- do.call(rbind, data)
        
        Heatmap(
            matrix=data, name="expression\nby sample",
            show_column_names=FALSE, cluster_columns=FALSE, 
            show_row_names=TRUE, row_names_side="right")
        
        # row annotation: adj. p-values
        cols <- colnames(top) %in% c("FDR", "adj.P.Val", "p_adj")
        s <- top[, cols] <= threshold
        s[is.na(s)] <- FALSE
        
        df_s <- data.frame(cluster_id=top$cluster_id, s=as.numeric(s)) 
        if (type == "DS") df_s$marker <- top$marker
        
        row_anno_df <- data.frame("significant"=factor(df_s$s, 
            levels=c(0,1), labels=c("no","yes")), check.names=FALSE)
        row_anno <- rowAnnotation(df=row_anno_df, width=unit(.8, "cm"),
            col=list("significant"=c("no"="gray90", "yes"="red2")))
        
        # combine panels
        main <- switch(type, DA="top DA clusters", 
            DS="top DS cluster-marker combinations")
        draw(hm1+hm2+hm3+row_anno, column_title=main, 
            column_title_gp=gpar(fontface="bold", fontsize=12))
    }
)