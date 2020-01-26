# ==============================================================================
# Heatmap for differental abundance & state analysis
# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
#' @title Plot differential heatmap
#' @description 
#' Heatmaps summarizing differental abundance 
#' & differential state testing results.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
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
#' @param col_anno
#'   logical. Should column annotations for each factor 
#'   in \code{metadata(x)} be included?
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
#' @author Lukas M Weber & 
#' Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
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
#' 
#' da <- diffcyt(sce, design = design, contrast = contrast, 
#'     analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
#'     clustering_to_use = "meta20")
#'     
#' ds <- diffcyt(sce, design = design, contrast = contrast, 
#'     analysis_type = "DS", method_DS = "diffcyt-DS-limma",
#'     clustering_to_use = "meta20")
#'     
#' # display test results for
#' # - top DA clusters
#' # - top DS cluster-marker combintations
#' plotDiffHeatmap(sce, da)
#' plotDiffHeatmap(sce, ds)
#' 
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom data.table data.table
#' @importFrom dplyr mutate_if
#' @importFrom methods is
#' @importFrom purrr map_depth
#' @importFrom scales hue_pal
#' @importFrom stats quantile
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors metadata
#' @export

plotDiffHeatmap <- function(x, y, 
    top_n = 20, all = FALSE, order = TRUE,
    th = 0.1, hm1 = TRUE, normalize = TRUE, 
    row_anno = TRUE, col_anno = TRUE) {
    
    .check_sce(x)
    es <- assay(x, "exprs")
    
    stopifnot(
        is.numeric(top_n), length(top_n) == 1,
        is.logical(order), length(order) == 1,
        is.numeric(th), length(th) == 1,
        is.logical(hm1), length(hm1) == 1,
        is.logical(normalize), length(normalize) == 1,
        is.logical(row_anno), length(row_anno) == 1)
    
    stopifnot(!is.null(k <- metadata(y$res)$clustering_name))
    k <- .check_validity_of_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
    factors <- select(
        as.data.frame(colData(x)), 
        -c("sample_id", "cluster_id"))
    
    y <- rowData(y$res)
    type <- .get_dt_type(y)
    
    # get clusters/cluster-marker combinations to plot
    if (order) y <- y[order(y$p_adj), , drop = FALSE]
    if (all | top_n > nrow(y)) top_n <- nrow(y)
    top <- as.data.frame(y[seq_len(top_n), ])
    top <- mutate_if(top, is.factor, as.character)
    
    # 1st heatmap: median type-marker expression by cluster
    if (hm1) {
        ms_by_k <- t(.agg(x[type_markers(x)], "cluster_id"))[top$cluster_id, ]
        qs <- quantile(ms_by_k, probs = c(.01, .5, .99), na.rm = TRUE)
        hm_cols <- colorRamp2(qs, c("royalblue3", "white", "tomato2"))
        hm1 <- .diff_hm(ms_by_k, hm_cols, "expression",
            cluster_rows = !order, xlab = "type_markers",
            row_title = "cluster_id"[!is.null(hm1)], 
            row_names_side = "left")
    } else {
        hm1 <- NULL
    }
    
    # column annotation: factors
    if (col_anno) {
        m <- match(levels(x$sample_id), x$sample_id)
        df <- data.frame(factors[m, ], row.names = NULL)
        col_anno <- .anno_factors(df, "column")
    } else {
        col_anno <- NULL
    }
    
    # 2nd heatmap:
    if (type == "DA") {
        # relative cluster abundances by sample
        cnts <- table(x$cluster_id, x$sample_id)
        frqs <- prop.table(cnts, 2)
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
            top_annotation = col_anno,
            xlab = "sample_id", 
            row_title = "cluster_id")
    } else {
        # median state-marker expression by sample
        cs_by_ks <- .split_cells(x, c("cluster_id", "sample_id"))
        ms_by_ks <- t(mapply(function(k, g)
            vapply(cs_by_ks[[k]], function(cs)
                median(es[g, cs, drop=FALSE]),
                numeric(1)),
            k = top$cluster_id, 
            g = top$marker_id))
        if (!is.null(hm1)) {
            rownames(ms_by_ks) <- top$marker_id 
        } else { 
            rownames(ms_by_ks) <- sprintf("%s(%s)", 
                top$marker_id, top$cluster_id)
        }
        if (normalize) 
            ms_by_ks <- .z_normalize(ms_by_ks) 
            hm2 <- .diff_hm(matrix=ms_by_ks, cluster_rows=!order,
                name=paste0("normalized\n"[normalize], "expression"),
                col=c("skyblue", "cornflowerblue", "royalblue", 
                    "black", "orange3", "orange", "gold"),
                xlab="sample_id", top_annotation=col_anno,
                row_names_side=c("right", "left")[as.numeric(is.null(hm1)) + 1])
    }
    
    # row annotation: significant = (adj. p-values <= th)
    if (row_anno) {
        s <- top$p_adj <= th
        s[is.na(s)] <- FALSE
        s <- as.matrix(c("no", "yes")[as.numeric(s)+1])
        rownames(s) <- format(top$p_adj, digits = 2)
        row_anno <- Heatmap(
            matrix = s, name = "significant",
            col = c(no="lightgrey", yes="limegreen"), 
            width = unit(5, "mm"), rect_gp = gpar(col = "white"),
            show_row_names = TRUE, row_names_side = "right")
    } else {
        row_anno <- NULL
    }
    
    # combine panels
    main <- switch(type, 
        DA = "top DA clusters", 
        DS = "top DS cluster-marker combinations")
    suppressWarnings(draw(hm1 + hm2 + row_anno, 
        column_title = main, auto_adjust = FALSE,
        column_title_gp = gpar(fontface = "bold", fontsize = 12)))
}