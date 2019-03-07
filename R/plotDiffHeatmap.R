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
#' @param row_anno
#'   logical. Should a row annotation indicating whether cluster (DA) 
#'   or cluster-marker combinations (DS) are significant, 
#'   as well as adjusted p-values be included?
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
#' @importFrom scales hue_pal
#' @importFrom stats quantile
# ------------------------------------------------------------------------------

setMethod(f = "plotDiffHeatmap",
    signature = signature(x = "matrix", y = "SummarizedExperiment"),
    definition = function(x, y, 
        top_n = 20, all = FALSE, order = TRUE,
        th = 0.1, hm1 = TRUE, normalize = TRUE, 
        row_anno = TRUE, ...) {
        
        # validity checks
        stopifnot(length(top_n) == 1, is.numeric(top_n))
        stopifnot(length(all) == 1, is.logical(all))
        stopifnot(length(order) == 1, is.logical(order))
        stopifnot(length(th) == 1, is.numeric(th))
        stopifnot(length(hm1) == 1, is.logical(hm1))
        stopifnot(length(normalize) == 1, is.logical(normalize))
        stopifnot(length(row_anno) == 1, is.logical(row_anno))
        
        z <- list(...)
        sample_ids <- z$sample_ids
        cluster_ids <- z$cluster_ids
        marker_classes <- z$marker_classes
        factors <- z$factors
        
        y <- rowData(y)
        analysis_type <- .get_dt_type(y)
        
        # get clusters/cluster-marker combinations to plot
        if (order)
            y <- y[order(y$p_adj), , drop = FALSE]
        if (all | top_n > nrow(y)) 
            top_n <- nrow(y)
        top <- y[seq_len(top_n), ]
        
        # 1st heatmap: median type-marker expression by cluster
        if (hm1) {
            type_markers <- colnames(x)[marker_classes == "type"]
            meds <- .calc_meds(x[, type_markers], "c", cluster_ids, NULL, top)
            qs <- quantile(meds, probs = c(.01, .5, .99), na.rm = TRUE)
            hm_cols <- colorRamp2(qs, c("royalblue3", "white", "tomato2"))
            hm1 <- .diff_hm(
                matrix = meds, 
                col = hm_cols, 
                name = "expression",
                cluster_rows = !order,
                xlab = "type_markers", 
                row_title = "cluster_id", 
                row_names_side = "left")
        } else {
            hm1 <- NULL
        }
        
        # column annotation: factors
        m <- match(levels(sample_ids), sample_ids)
        df <- data.frame(factors[m, ], row.names = NULL)
        lvls <- lapply(seq_len(ncol(df)), function(i) levels(df[[i]]))
        nlvls <- vapply(lvls, length, numeric(1))
        cols <- hue_pal()(sum(nlvls))
        names(cols) <- unlist(lvls)
        cols <- split(cols, rep.int(seq_len(ncol(df)), nlvls))
        names(cols) <- names(df)
        col_anno <- columnAnnotation(df, col = cols, gp = gpar(col = "white"))
        
        # 2nd heatmap:
        if (analysis_type == "DA") {
            # relative cluster abundances by sample
            frqs <- prop.table(table(cluster_ids, sample_ids), 2)
            frqs <- frqs[top$cluster_id, ]
            frqs <- as.matrix(unclass(frqs))
            if (normalize) {
                frqs <- .z_normalize(asin(sqrt(frqs)))
                at <- seq(-2.5, 2.5, 0.5)
                labels <- at
                labels[-seq(2, length(at), 2)] <- ""
            } else {
                min <- floor(min(frqs)/0.1)*0.1
                max <- ceiling(max(frqs)/0.1)*0.1
                at <- seq(min, max, 0.1)
                labels <- at
            }
            hm2 <- .diff_hm(matrix = frqs, cluster_rows = !order, 
                col = c("skyblue", "cornflowerblue", "royalblue", 
                    "black", "orange3", "orange", "gold"),
                name = paste0("normalized\n"[normalize], "frequency"),
                show_row_names = is.null(hm1), row_names_side = "left",
                heatmap_legend_param = list(at = at, labels = labels),
                xlab = "sample_id", top_annotation = col_anno)
        } else {
            # median state-marker expression by sample
            meds <- .calc_meds(x, "cs", cluster_ids, sample_ids, top)
            if (normalize) meds <- .z_normalize(meds) 
            hm2 <- .diff_hm(matrix = meds, cluster_rows=!order,
                name = paste0("normalized\n"[normalize], "expression"),
                col = c("skyblue", "cornflowerblue", "royalblue", 
                    "black", "orange3", "orange", "gold"),
                xlab = "sample_id", top_annotation = col_anno)
        }
        
        # row annotation: significant = (adj. p-values <= th)
        if (row_anno) {
            s <- top$p_adj <= th
            s[is.na(s)] <- FALSE
            s <- as.matrix(c("no", "yes")[as.numeric(s)+1])
            rownames(s) <- format(top$p_adj, digits = 2)
            row_anno <- Heatmap(
                matrix = s, 
                name = "significant",
                col = c(no = "lightgrey", yes = "limegreen"),
                width = unit(5, "mm"),
                rect_gp = gpar(col = "white"),
                show_row_names = TRUE,
                row_names_side = "right")
        } else {
            row_anno <- NULL
        }
        
        # combine panels
        main <- switch(analysis_type, 
            DA = "top DA clusters", 
            DS = "top DS cluster-marker combinations")
        suppressWarnings(
            draw(hm1 + hm2 + row_anno, column_title = main,
                column_title_gp = gpar(fontface = "bold", fontsize = 12)))
    }
)

# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
#' @importFrom dplyr %>% select
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom S4Vectors metadata
setMethod(f="plotDiffHeatmap",
    signature=signature(x="daFrame", y="SummarizedExperiment"),
    definition=function(x, y, top_n = 20, all = FALSE, order = TRUE,
        th = 0.1, hm1 = TRUE, normalize = TRUE, row_anno = TRUE, ...) {
        
        # get cluster IDs
        k <- metadata(y)$clustering_name
        k <- .check_validity_of_k(x, k)
        cluster_ids <- .get_cluster_ids(x, k)

        plotDiffHeatmap(exprs(x), y,
            top_n, all, order, th, hm1, normalize, row_anno,
            sample_ids=sample_ids(x), 
            cluster_ids=cluster_ids,
            marker_classes=colData(x)$marker_class,
            factors=rowData(x) %>% data.frame %>% 
                select(-c("sample_id", "cluster_id")))
    }
)

# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
#' @importFrom dplyr %>% select
#' @importFrom SummarizedExperiment assay colData rowData
setMethod(f="plotDiffHeatmap",
    signature=signature(x="SummarizedExperiment", y="SummarizedExperiment"),
    definition=function(x, y, top_n = 20, all = FALSE, order = TRUE,
        th = 0.1, hm1 = TRUE, normalize = TRUE, row_anno = TRUE, ...) {
        
        plotDiffHeatmap(assay(x), y, 
            top_n, all, order, th, hm1, normalize, row_anno,
            sample_ids=rowData(x)$sample_id,
            cluster_ids=rowData(x)$cluster_id,
            marker_classes=colData(x)$marker_class,
            factors=rowData(x) %>% data.frame %>% 
                select(-c("sample_id", "cluster_id")))
    }
)

# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
setMethod(f="plotDiffHeatmap", 
    signature=signature(x="ANY", y="list"), 
    definition=function(x, y, top_n = 20, all = FALSE, order = TRUE,
        th = 0.1, hm1 = TRUE, normalize = TRUE, row_anno = TRUE, ...) {
        
        if (all(c("res", "d_counts", "d_medians") %in% names(y))) {
            plotDiffHeatmap(x, y$res, 
                top_n, all, order, th, hm1, normalize, row_anno)
        } else {
            stop(deparse(substitute(y)), " does not seem to be ", 
                "a valid differential test result.\n",
                "Should be a 'SummarizedExperiment' as returned by ", 
                "'diffcyt::testDA_*()' or 'diffcyt::testDS_*()'.")
        }
    }
)