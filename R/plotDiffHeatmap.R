# ==============================================================================
# Heatmap for differental abundance & state analysis
# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
#' @title Plot differential heatmap
#' @description 
#' Heatmaps summarizing differental abundance 
#' & differential state testing results.
#' 
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param y 
#'   a \code{SummarizedExperiment} containing differential testing
#'   results as returned by one of \code{\link[diffcyt]{testDA_edgeR}}, 
#'   \code{\link[diffcyt]{testDA_voom}}, \code{\link[diffcyt]{testDA_GLMM}}, 
#'   \code{\link[diffcyt]{testDS_limma}}, or \code{\link[diffcyt]{testDS_LMM}}.
#' @param top_n 
#'   numeric. Number of top clusters (if \code{type = "DA"}) or
#'   cluster-marker combinations (if \code{type = "DS"}) to display.
#' @param all 
#'   logical. Specifies whether all clusters or cluster-marker combinations 
#'   should be displayed. If \code{TRUE}, \code{top_n} will be ignored.
#' @param order 
#'   logical. Should results be ordered by significance?
#' @param th 
#'   numeric. Threshold on adjusted p-values below which clusters (DA) 
#'   or cluster-marker combinations (DS) should be considered significant.
#' 
#' @details 
#' For DA tests, \code{plotDiffHeatmap} will display
#' \itemize{
#' \item median (arcsinh-transformed) cell-type marker expressions (across all samples)
#' \item cluster abundances by samples
#' \item row annotations indicating the significance of detected clusters
#' }
#' For DS tests, \code{plotDiffHeatmap} will display
#' \itemize{
#' \item median (arcsinh-transformed) cell-type marker expressions (across all samples)
#' \item median (arcsinh-transformed) cell-state marker expressions (across all samples)
#' \item median (arcsinh-transformed) cell-state marker expressions by sample
#' \item row annotations indicating the significance of detected cluster-marker combinations
#' }
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @author Lukas M Weber, Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @import ComplexHeatmap 
#' @importFrom circlize colorRamp2
#' @importFrom dplyr group_by_ summarize_at
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @importFrom tidyr complete
# ------------------------------------------------------------------------------

setMethod(f="plotDiffHeatmap", 
    signature=signature(x="daFrame", y="SummarizedExperiment"), 
    definition=function(x, y, top_n=20, all=FALSE, order=TRUE, th=0.1) {

        # get differential analysis type
        type <- get_dt_type(y)
        
        # get no. of clusters
        k <- switch(type, 
            DA = nrow(y),
            DS = nrow(y) / nlevels(rowData(y)$marker))
        
        # validity check
        check_daFrame_dt <- function(x, y) {
            
        }
        
        # get cluster IDs
        cluster_ids <- cluster_codes(x)[, k][cluster_ids(x)]
        df <- data.frame(exprs(x), 
            sample_id=sample_ids(x), 
            cluster_id=cluster_ids)
        
        # compute medians across samples & clusters
        meds_by_sample <- df %>% group_by_(~sample_id) %>% 
            summarize_at(colnames(x), funs(median))
        meds_by_cluster <- df %>% group_by_(~cluster_id) %>% 
            summarise_at(colnames(x), funs(median))

        # color scale: 
        # 1%, 50%, 99% percentiles across medians & markers
        qs <- quantile(unlist(meds_by_cluster), c(.01, .5, .99), TRUE)
        hm_cols <- colorRamp2(qs, c("royalblue3", "white", "tomato2"))
        
        # get clusters or cluster-marker combinations to plot
        y <- rowData(y)
        if (order) {
            cols <- colnames(y) %in% c("FDR", "adj.P.Val", "p_adj")
            y <- y[order(y[, cols]), , drop=FALSE]
        }
        if (all | top_n > nrow(y)) 
            top_n <- nrow(y)
        top <- y[seq_len(top_n), ]
        meds_by_cluster <- meds_by_cluster[
            match(top$cluster_id, rownames(meds_by_cluster)), ]
        
        # 1st heatmap:
        # median cell-type marker expressions across clusters
        hm1 <- Heatmap(
            meds_by_cluster[, type_markers(x)], hm_cols, "expression", 
            column_title="type_markers", column_title_side="bottom", 
            cluster_columns=FALSE, row_names_side="left", 
            clustering_distance_rows="euclidean",
            clustering_method_rows="median")
        
        meds <- df %>% 
            group_by_(~cluster_id, ~sample_id) %>% 
            summarise_all(funs(median))
        meds <- lapply(colnames(x), function(marker) 
            acast(meds, cluster_id~sample_id, value.var=marker, fill=0))
        meds <- setNames(meds, colnames(x))
        meds <- meds[top$cluster_id, , drop=FALSE]
       
        # 2nd heatmap:
        # type = "DA": cluster sizes across samples
        # type = "DS": median cell-state marker expressions across clusters
        if (type == "DA") {
            counts <- df %>% count(cluster_id, sample_id) %>% complete(sample_id)
            counts <- reshape2::acast(counts, cluster_id~sample_id, value.var="n", fill=0)
            counts <- counts[top$cluster_id, , drop=FALSE]
            
            hm2 <- Heatmap(
                counts, colorRamp2(range(counts), c("black", "gold")), 
                "n_cells", cluster_columns=FALSE, show_row_names=FALSE)
        } else if (type == "DS") { 
            hm2 <- Heatmap(
                meds_by_cluster[, state_markers(x)], hm_cols, 
                column_title="state_markers", column_title_side="bottom", 
                cluster_columns=FALSE, row_names_side = "left", 
                clustering_distance_rows="euclidean",
                clustering_method_rows="median",
                show_heatmap_legend=FALSE)
        }
        
        # if type = "DS", 3rd heatmap:
        # median cell-state marker expressions across samples
        
        # subset top markers & clusters
        data <- mapply(function(marker, ids) marker[ids, , drop = FALSE], 
            meds[top$marker], top$cluster_id, SIMPLIFY=FALSE)
        data <- do.call(rbind, data)
        
        hm3 <- Heatmap(
            matrix=data, name="expression\nby sample",
            show_column_names=FALSE, cluster_columns=FALSE, 
            show_row_names=TRUE, row_names_side="right")
        
        # row annotation: adj. p-values
        cols <- colnames(top) %in% c("FDR", "adj.P.Val", "p_adj")
        s <- top[, cols] <= th
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

# ==============================================================================
# method for when 'y' is a list as returned by diffcyt()
# ------------------------------------------------------------------------------
setMethod(f="plotDiffHeatmap", 
    signature=signature(x="daFrame", y="list"), 
    definition=function(x, y, top_n=20, all=FALSE, order=TRUE, th=0.1) {
        if (all.equal(names(y), c("res", "d_se", "d_counts", "d_medians", 
            "d_medians_by_cluster_marker", "d_medians_by_sample_marker"))) {
            plotDiffHeatmap(x, y$res, top_n=20, all=FALSE, order=TRUE, th=0.1)
        } else {
            stop(deparse(substitute(y)), " does not seem to be ", 
                "a valid differential test result.\n",
                "Should be a 'SummarizedExperiment' as returned by ", 
                "'diffcyt::testDA_*()' or 'diffcyt::testDS_*()'.")
        }
    }
)