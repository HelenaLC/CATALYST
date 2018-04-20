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
#' \item median (arcsinh-transformed) cell-state marker expressions by sample
#' \item row annotations indicating the significance of detected cluster-marker combinations
#' }
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @author Lukas M Weber, Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @examples
#' 
#' 
#' @import ComplexHeatmap 
#' @importFrom circlize colorRamp2
#' @importFrom dplyr group_by_ summarize_all count
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @importFrom tidyr complete
#' @importFrom reshape2 acast
# ------------------------------------------------------------------------------

setMethod(f="plotDiffHeatmap", 
    signature=signature(x="daFrame", y="SummarizedExperiment"), 
    definition=function(x, y, top_n=20, all=FALSE, order=TRUE, th=0.1) {

        # get differential analysis type
        y <- rowData(y)
        type <- get_dt_type(y)
        
        # get no. of clusters & cluster IDs
        k <- switch(type, 
            DA = nrow(y),
            DS = nrow(y) / nlevels(y$marker))
        k <- check_validity_of_k(x, k)
        cluster_ids <- cluster_codes(x)[, k][cluster_ids(x)]
        
        # compute medians by samples & clusters
        df <- data.frame(exprs(x), sample_id=sample_ids(x), cluster_id=cluster_ids)
        meds_by_sample <- data.frame(df %>% group_by_(~sample_id) %>% 
                summarize_at(colnames(x), median), row.names=1)
        meds_by_cluster <- data.frame(df %>% group_by_(~cluster_id) %>% 
                summarise_at(colnames(x), median), row.names=1)

        # color scale: 1%, 50%, 99% percentiles
        qs <- quantile(meds_by_cluster, c(.01, .5, .99), TRUE)
        hm_cols <- colorRamp2(qs, c("royalblue3", "white", "tomato2"))
        
        # get clusters or cluster-marker combinations to plot
        if (order)
            y <- y[order(y$p_adj), , drop=FALSE]
        if (all | top_n > nrow(y)) 
            top_n <- nrow(y)
        top <- y[seq_len(top_n), ]
        meds_by_cluster <- meds_by_cluster[top$cluster_id, , drop=FALSE]
        meds_by_cluster <- as.matrix(meds_by_cluster)
        rownames(meds_by_cluster) <- top$cluster_id
        
        # 1st heatmap:
        # median type-marker expressions by cluster
        hm1 <- diff_hm(matrix=meds_by_cluster[, type_markers(x)], 
            col=hm_cols, name="expression", xlab="type_markers", 
            row_title="cluster_id", row_names_side="left")
        
        # 2nd heatmap:
        if (type == "DA") {
            # cluster sizes by sample
            n_cells <- df %>% count(cluster_id, sample_id) %>% complete(sample_id)
            n_cells <- reshape2::acast(n_cells, cluster_id~sample_id, value.var="n", fill=0)
            n_cells <- n_cells[top$cluster_id, , drop=FALSE]
            
            hm2 <- diff_hm(matrix=n_cells, name="n_cells",
                col=colorRamp2(range(n_cells), c("navy", "yellow")),
                xlab="sample_id", show_row_names=FALSE)
        } else if (type == "DS") { 
            # median state-marker expression by sample
            meds <- df %>% 
                group_by_(~cluster_id, ~sample_id) %>% 
                summarise_all(funs(median))
            meds <- lapply(colnames(x), function(marker) 
                acast(meds, cluster_id~sample_id, value.var=marker, fill=0))
            meds <- setNames(meds, colnames(x))
            meds <- mapply(function(marker, ids) marker[ids, , drop = FALSE], 
                meds[top$marker], top$cluster_id, SIMPLIFY=FALSE)
            meds <- do.call(rbind, meds)
            rownames(meds) <- paste0(top$marker, sprintf("(%s)", top$cluster_id))

            hm2 <- diff_hm(matrix=meds, name="expression\nby sample",
                col=colorRamp2(range(meds, na.rm=TRUE), c("navy", "yellow")),
                xlab="sample_id")
        }

        # row annotation: significant = (adj. p-values <= th)
        s <- top[, "p_adj"] <= th
        s[is.na(s)] <- FALSE
        df_s <- data.frame(cluster_id=top$cluster_id, s=as.numeric(s)) 
        if (type == "DS") df_s$marker <- top$marker
        row_anno_df <- data.frame("significant"=factor(df_s$s, 
            levels=c(0,1), labels=c("no","yes")), check.names=FALSE)
        
        row_anno <- rowAnnotation(df=row_anno_df, 
            gp=gpar(col="white"), width=unit(.4, "cm"),
            col=list("significant"=c("no"="gray90", "yes"="limegreen")))
        
        # combine panels
        main <- switch(type, 
            DA = "top DA clusters", 
            DS = "top DS cluster-marker combinations")
        draw(hm1 + hm2 + row_anno, column_title=main,
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