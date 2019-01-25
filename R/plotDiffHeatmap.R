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
#'   a \code{\link{daFrame}} or \code{SummarizedExperiment}.
#' @param y 
#'   a \code{SummarizedExperiment} containing differential testing
#'   results as returned by one of \code{\link[diffcyt]{testDA_edgeR}}, 
#'   \code{\link[diffcyt]{testDA_voom}}, \code{\link[diffcyt]{testDA_GLMM}}, 
#'   \code{\link[diffcyt]{testDS_limma}}, or \code{\link[diffcyt]{testDS_LMM}}.
#'   Alternatively, a list as returned by \code{\link[diffcyt]{diffcyt}}.
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
#' @param hm1 
#'   logical. Specifies whether the left-hand side heatmap should be plotted.
#' @param normalize
#'   logical. Specifies whether Z-score normalized values should be plotted 
#'   in the right-hand side heatmap. If \code{y} contains DA analysis results, 
#'   relative population abundances will be arcsine-square-root scaled 
#'   prior to normalization.
#' @param anno_signif
#'   logical. Should a row annotation indicating whether cluster (DA) 
#'   or cluster-marker combinations (DS) are significant be included?
#' @param anno_padj
#'   logical. Should a row annotation for the adjusted p-values of 
#'   each cluster (DA) or cluster-marker combinations (DS) be included?
#' 
#' @details 
#' For DA tests, \code{plotDiffHeatmap} will display
#' \itemize{
#'   \item{median (arcsinh-transformed) 
#'     cell-type marker expressions (across all samples)}
#'   \item{cluster abundances by samples}
#'   \item{row annotations indicating if detected clusteres
#'     are significant (i.e. adj. p-value >= \code{th})}
#' }
#' For DS tests, \code{plotDiffHeatmap} will display
#'   \itemize{
#'   \item{median (arcsinh-transformed) 
#'     cell-type marker expressions (across all samples)}
#'   \item{median (arcsinh-transformed) 
#'     cell-state marker expressions by sample}
#'   \item{row annotations indicating if detected cluster-marker combinations
#'     are significant (i.e. adj. p-value >= \code{th})}
#' }
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @author Lukas M Weber and 
#' Helena Lucia Crowell \email{helena.crowell@uzh.ch}
#' 
#' @examples
#' # construct daFrame
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run clustering
#' re <- cluster(re)
#' 
#' ## differential analysis
#' library(diffcyt)
#' 
#' # create design & constrast matrix
#' design <- createDesignMatrix(PBMC_md, cols_design=3:4)
#' contrast <- createContrast(c(0, 1, 0, 0, 0))
#' 
#' # test for
#' # - differential abundance (DA) of clusters
#' # - differential states (DS) within clusters
#' da <- diffcyt(re, design = design, contrast = contrast, 
#'     analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
#'     clustering_to_use = "meta20")
#' ds <- diffcyt(re, design = design, contrast = contrast, 
#'     analysis_type = "DS", method_DS = "diffcyt-DS-limma",
#'     clustering_to_use = "meta20")
#'     
#' # display test results for
#' # - top DA clusters
#' # - top DS cluster-marker combintations
#' plotDiffHeatmap(re, da)
#' plotDiffHeatmap(re, ds)
#' 
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom dplyr group_by_ count summarise_all summarise_at 
#' @importFrom magrittr %>%
#' @importFrom Matrix colSums
#' @importFrom stats quantile
#' @importFrom tidyr complete
#' @importFrom reshape2 acast
# ------------------------------------------------------------------------------

setMethod(f="plotDiffHeatmap", 
    signature=signature(x="matrix", y="SummarizedExperiment"), 
    definition=function(x, y, top_n=20, all=FALSE, 
        order=TRUE, th=0.1, hm1=TRUE, normalize=FALSE, 
        anno_signif=TRUE, anno_padj=TRUE, ...) {
        
        # validity checks
        stopifnot(length(top_n) == 1, is.numeric(top_n))
        stopifnot(length(all) == 1, is.logical(all))
        stopifnot(length(order) == 1, is.logical(order))
        stopifnot(length(th) == 1, is.numeric(th))
        stopifnot(length(hm1) == 1, is.logical(hm1))
        stopifnot(length(normalize) == 1, is.logical(normalize))
        
        z <- list(...)
        sample_ids <- z$sample_ids
        cluster_ids <- z$cluster_ids
        marker_classes <- z$marker_classes
        
        # get differential analysis type
        y <- rowData(y)
        type <- get_dt_type(y)
        
        # get clusters or cluster-marker combinations to plot
        if (order)
            y <- y[order(y$p_adj), , drop=FALSE]
        if (all | top_n > nrow(y)) 
            top_n <- nrow(y)
        top <- y[seq_len(top_n), ]
        
        df <- data.frame(x, sample_id=sample_ids, cluster_id=cluster_ids)
        if (hm1) {
            meds_by_cluster <- data.frame(df %>% group_by_(~cluster_id) %>% 
                    summarise_at(colnames(x), median), row.names=1)
            if (nlevels(cluster_ids) == 1)
                meds_by_cluster <- meds_by_cluster[, -1]
            
            # color scale: 1%, 50%, 99% percentiles
            qs <- quantile(meds_by_cluster[, marker_classes != "none"], c(.01, .5, .99), TRUE)
            hm_cols <- circlize::colorRamp2(qs, c("royalblue3", "white", "tomato2"))
            
            meds_by_cluster <- meds_by_cluster[top$cluster_id, , drop=FALSE]
            meds_by_cluster <- as.matrix(meds_by_cluster)
            rownames(meds_by_cluster) <- top$cluster_id
            
            # 1st heatmap:
            # median type-marker expressions by cluster
            hm1 <- diff_hm(cluster_rows=!order,
                matrix=meds_by_cluster[, marker_classes == "type"], 
                col=hm_cols, name="expression", xlab="type_markers", 
                row_title="cluster_id", row_names_side="left")
        } else {
            hm1 <- NULL
        }
        
        # 2nd heatmap:
        if (type == "DA") {
            # compute relative cluster sizes by sample
            freqs <- prop.table(table(df$cluster_id, df$sample_id), 2)
            freqs <- freqs[top$cluster_id, , drop=FALSE]
            freqs <- as.matrix(unclass(freqs))
            if (normalize)
                freqs <- z_normalize(asin(sqrt(freqs)))
            hm2 <- diff_hm(matrix=freqs, cluster_rows=!order,
                name=paste0("normalized\n"[normalize], "frequency"),
                col=colorRamp2(range(freqs), c("navy", "yellow")),
                xlab="sample_id", show_row_names=FALSE)
        } else if (type == "DS") { 
            # median state-marker expression by sample
            meds <- df %>% 
                group_by_(~cluster_id, ~sample_id) %>% 
                summarise_all(funs(median))
            meds <- lapply(colnames(x), function(marker_id) 
                acast(meds, cluster_id~sample_id, 
                    value.var=marker_id, fill=0))
            meds <- setNames(meds, colnames(x))
            meds <- mapply(function(marker, ids) marker[ids, , drop = FALSE], 
                meds[top$marker_id], top$cluster_id, SIMPLIFY=FALSE)
            meds <- do.call(rbind, meds)
            if (normalize)
                meds <- z_normalize(meds) 
            rownames(meds) <- paste0(top$marker, sprintf("(%s)", top$cluster_id))
            hm2 <- diff_hm(matrix=meds, cluster_rows=!order,
                name=paste0("normalized\n"[normalize], "expression"),
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
        
        if (anno_signif) {
            row_anno <- rowAnnotation(df=row_anno_df, 
                gp=gpar(col="white"), width=unit(.4, "cm"),
                col=list("significant"=c("no"="gray90", "yes"="limegreen")))
        } else {
            row_anno <- NULL
        }
        
        if (anno_padj) {
            pval_hm <- Heatmap(name="p_adj",
                matrix=matrix(top$p_adj, ncol=1, dimnames=list(
                    format(top$p_adj, digits=3, scientific=TRUE), NULL)), 
                col=colorRamp2(range(top$p_adj, na.rm=TRUE), c("black", "lightgrey")),
                width=unit(.4, "cm"), rect_gp=gpar(col="white"), 
                show_row_names=TRUE, row_names_side="right", 
                cluster_rows=FALSE, show_column_names=FALSE)
        } else {
            pval_hm <- NULL
        }
        
        # combine panels
        main <- switch(type, 
            DA = "top DA clusters", 
            DS = "top DS cluster-marker combinations")
        suppressWarnings(
            draw(hm1 + hm2 + row_anno + pval_hm, column_title=main,
                column_title_gp=gpar(fontface="bold", fontsize=12)))
    }
)

# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
setMethod(f="plotDiffHeatmap",
    signature=signature(x="daFrame", y="SummarizedExperiment"),
    definition=function(x, y, top_n=20, all=FALSE, 
        order=TRUE, th=0.1, hm1=TRUE, normalize=FALSE, ...) {
        
        # get cluster IDs
        k <- metadata(y)$clustering_name
        k <- check_validity_of_k(x, k)
        cluster_ids <- cluster_codes(x)[, k][cluster_ids(x)]
        
        plotDiffHeatmap(exprs(x), y, top_n, all, order, th, hm1, normalize,
            sample_ids=sample_ids(x), 
            cluster_ids=cluster_ids,
            marker_classes=colData(x)$marker_class)
    }
)

# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
setMethod(f="plotDiffHeatmap",
    signature=signature(x="SummarizedExperiment", y="SummarizedExperiment"),
    definition=function(x, y, top_n=20, all=FALSE, 
        order=TRUE, th=0.1, hm1=TRUE, normalize=FALSE, ...) {
        
        plotDiffHeatmap(assay(x), y, top_n, all, order, th, hm1, normalize,
            sample_ids=rowData(x)$sample_id,
            cluster_ids=rowData(x)$cluster_id,
            marker_classes=colData(x)$marker_class)
    }
)

# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
setMethod(f="plotDiffHeatmap", 
    signature=signature(x="ANY", y="list"), 
    definition=function(x, y, top_n=20, all=FALSE, 
        order=TRUE, th=0.1, hm1=TRUE, normalize=FALSE, ...) {
        
        if (all(c("res", "d_counts", "d_medians") %in% names(y))) {
            plotDiffHeatmap(x, y$res, top_n, all, order, th, hm1, normalize)
        } else {
            stop(deparse(substitute(y)), " does not seem to be ", 
                "a valid differential test result.\n",
                "Should be a 'SummarizedExperiment' as returned by ", 
                "'diffcyt::testDA_*()' or 'diffcyt::testDS_*()'.")
        }
    }
)