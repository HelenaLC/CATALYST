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
#' @param assay character string specifying which assay 
#'   data to use; valid values are \code{assayNames(x)}.
#' @param fun character string specifying the function to use 
#'   as summary statistic for aggregation of expression values.
#' @param scale logical specifying whether expression values should be scaled
#'   between 0 and 1 using lower (1\%) and upper (99\%) quantiles as boundaries.
#' @param normalize
#'   logical. Specifies whether Z-score normalized values should be plotted 
#'   in the right-hand side heatmap. If \code{y} contains DA analysis results, 
#'   relative population abundances will be arcsine-square-root scaled 
#'   prior to normalization.
#' @param row_anno logical specifying whether to invlude a row annotation 
#'   indicating whether cluster (DA) or cluster-marker combinations (DS) 
#'   are significant, labeled with adjusted p-values.
#' @param col_anno logical specifying whether to include column annotations 
#'   for all non-numeric cell metadata variables; or a character vector 
#'   in \code{names(colData(x))} to include only a subset of annotations.
#'   (Only variables that map uniquely to each sample will be included)
#' @param hm_pal character vector of colors to interpolate for the heatmap. 
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
#' # - top DS cluster-marker combinations
#' plotDiffHeatmap(sce, da)
#' plotDiffHeatmap(sce, ds)
#' 
#' # visualize results for subset of clusters
#' sub <- filterSCE(sce, k = "meta20", cluster_id %in% seq_len(5))
#' plotDiffHeatmap(sub, da, order = FALSE)
#' 
#' # visualize results for selected feature
#' # & include only selected annotation
#' plotDiffHeatmap(sce["pS6", ], ds, 
#'   col_anno = "condition", all = TRUE)
#' 
#' @importFrom ComplexHeatmap rowAnnotation anno_simple row_anno_text Heatmap
#' @importFrom dplyr mutate_at mutate_if
#' @importFrom grid gpar unit
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales scientific
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom S4Vectors metadata
#' @export

plotDiffHeatmap <- function(x, y, 
    top_n = 20, all = FALSE, order = TRUE, th = 0.1, 
    assay = "exprs", fun = c("median", "mean", "sum"), 
    scale = TRUE, normalize = TRUE,
    col_anno = TRUE, row_anno = TRUE,
    hm_pal = rev(brewer.pal(11, ifelse(
        is.null(assayNames(y$res)), "RdYlBu", "RdBu")))) {
    
    # check validity of input arguments
    fun <- match.arg(fun)
    args <- as.list(environment())
    .check_args_plotDiffHeatmap(args)
    
    stopifnot(!is.null(k <- metadata(y$res)$clustering_name))
    k <- .check_k(x, k); m <- NULL
    x$cluster_id <- cluster_ids(x, k)
    
    y <- rowData(y$res)
    type <- .get_dt_type(y)
    
    # subset results in case input SCE has been filtered
    i <- y$cluster_id %in% levels(x$cluster_id)
    if (type == "DS") i <- i & y$marker_id %in% rownames(x)
    y <- y[i, , drop = FALSE]

    if (nrow(y) == 0)
        stop("No results remaining;",
            " perhaps 'x' has been filtered?")
    
    # get clusters/cluster-marker combinations to plot
    if (order) y <- y[order(y$p_adj), , drop = FALSE]
    if (all || top_n > nrow(y)) top_n <- nrow(y)
    top <- data.frame(y[seq_len(top_n), ])
    top <- mutate_if(top, is.factor, as.character)
    
    # column annotation of non-numeric cell metadata variables
    if (!isFALSE(col_anno)) {
        top_anno <- .anno_factors(x, levels(x$sample_id), col_anno, "column")
    } else top_anno <- NULL
    
    # row annotation: significant = (adj. p-values <= th)
    if (row_anno) {
        s <- top$p_adj <= th
        s[is.na(s)] <- FALSE
        ss <- c("no", "yes")
        s <- factor(ss[s+1], ss)
        txt <- scientific(top$p_adj, 2)
        right_anno <- rowAnnotation(
            significant = anno_simple(s, 
                gp = gpar(col = "white"),
                simple_anno_size = unit(2, "mm"),
                col = c(no = "lightgrey", yes = "lightgreen")),
            "foo" = row_anno_text(txt,
                width = unit(4, "mm"),
                gp = gpar(fontsize = 8)),
            show_annotation_name = FALSE)
    } else right_anno <- NULL

    switch(type, 
        # relative cluster abundances by sample
        DA = {
            ns <- table(x$cluster_id, x$sample_id)
            fq <- prop.table(ns, 2)
            fq <- fq[top$cluster_id, ]
            y <- as.matrix(unclass(fq))
            if (normalize) y <- .z_normalize(asin(sqrt(y)))
            Heatmap(
                matrix = y, 
                name = paste0("normalized\n"[normalize], "frequency"),
                col = hm_pal,
                na_col = "lightgrey", 
                rect_gp = gpar(col = "white"),  
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                row_names_side = "left",
                top_annotation = top_anno,
                right_annotation = right_anno)
        },
        # median state-marker expression by sample
        DS = {
            y <- assay(x, assay)
            cs <- .split_cells(x, c("cluster_id", "sample_id"))
            z <- t(mapply(function(k, g)
                vapply(cs[[k]], function(cs) {
                    if (length(cs) == 0) return(NA)
                    get(fun)(y[g, cs, drop = FALSE])
                }, numeric(1)),
                k = top$cluster_id, 
                g = top$marker_id))
            rownames(z) <- sprintf("%s(%s)", top$marker_id, top$cluster_id)
            if (normalize) z <- .z_normalize(z) 
            Heatmap(
                matrix = z,
                name = paste0("z-normalized\n"[normalize], "expression"),
                col = hm_pal,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                top_annotation = top_anno,
                row_names_side = "left",
                rect_gp = gpar(col = "white"),
                right_annotation = right_anno)
    }) 
}
