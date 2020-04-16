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
#' @param fun character string specifying the function to use 
#'   as summary statistic for aggregation of expression values.
#' @param scale logical specifying whether expression values should be scaled
#'   between 0 and 1 using lower (1\%) and upper (99\%) quantiles as boundaries.
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
#' @param hm1_pal,hm2_pal character vector of 
#'   colors to interpolate for each heatmap(s).  
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
#' # visualize results for subset of clusters
#' sub <- filterSCE(sce, k = "meta20", cluster_id %in% seq_len(5))
#' plotDiffHeatmap(sub, da)
#' 
#' @import ComplexHeatmap
#' @importFrom data.table data.table
#' @importFrom dplyr mutate_if
#' @importFrom grid unit.c
#' @importFrom methods is
#' @importFrom purrr map_depth
#' @importFrom stats quantile
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors metadata
#' @export

plotDiffHeatmap <- function(x, y, 
    top_n = 20, all = FALSE, order = TRUE, th = 0.1, 
    hm1 = TRUE, fun = c("median", "mean"), 
    scale = TRUE, normalize = TRUE, 
    row_anno = TRUE, col_anno = TRUE,
    hm1_pal = rev(brewer.pal(11, "RdBu")),
    hm2_pal = rev(brewer.pal(11, "PuOr"))) {
    
    # check validity of input arguments
    fun <- match.arg(fun)
    args <- as.list(environment())
    .check_args_plotDiffHeatmap(args)
    
    stopifnot(!is.null(k <- metadata(y$res)$clustering_name))
    k <- .check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
    
    y <- rowData(y$res)
    type <- .get_dt_type(y)
    
    # subset results in case input SCE has been filtered
    y <- y[y$cluster_id %in% levels(x$cluster_id), , drop = FALSE]

    # get clusters/cluster-marker combinations to plot
    if (order) y <- y[order(y$p_adj), , drop = FALSE]
    if (all || top_n > nrow(y)) top_n <- nrow(y)
    top <- data.frame(y[seq_len(top_n), ])
    top <- mutate_if(top, is.factor, as.character)
    
    # 1st heatmap: median type-marker expression by cluster
    if (hm1) {
        es <- assay(x, "exprs")
        if (scale) es <- .scale_exprs(es, 1)
        z <- x; assay(z, "exprs") <- es
        z <- z[type_markers(x), ]
        z <- .agg(z, "cluster_id", fun)
        z <- t(z)[top$cluster_id, ]
        if (type == "DS" && row_anno)
            rownames(z) <- sprintf("%s(%s)", 
                top$marker_id, top$cluster_id)
        hm1 <- Heatmap(
            matrix = z,
            name = paste0("scaled\n"[scale], "expression"),
            col = hm1_pal,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            row_title = "cluster_id",
            column_title = "type_markers",
            row_names_side = "left",
            column_title_side = "bottom",
            clustering_distance_rows = "euclidean",
            clustering_method_rows = "median",
            #row_names_gp = gpar(fontsize = 8),
            #column_names_gp = gpar(fontsize = 8),
            rect_gp = gpar(col = "white"))
    } else hm1 <- NULL
    
    # column annotation: factors
    cd <- data.frame(colData(x))
    cd <- select_if(cd, is.factor)
    cols_keep <- setdiff(colnames(cd), 
        c("cluster_id", "sample_id"))
    if (col_anno && length(cols_keep) > 0) {
        m <- match(levels(x$sample_id), x$sample_id)
        df <- data.frame(cd[m, cols_keep, drop = FALSE])
        col_anno <- .anno_factors(df, "column")
    } else col_anno <- NULL
    
    # row annotation: significant = (adj. p-values <= th)
    if (row_anno) {
        s <- top$p_adj <= th
        s[is.na(s)] <- FALSE
        ss <- c("no", "yes")
        s <- factor(ss[s+1], ss)
        txt <- format(top$p_adj, digits = 2)
        right_anno <- rowAnnotation(
            df = data.frame(significant = s),
            col = list(significant = c(no = "lightgrey", yes = "lightgreen")),
            "foo" = row_anno_text(txt),#, gp = gpar(fontsize = 6)),
            gp = gpar(col = "white"),
            show_annotation_name = FALSE,
            annotation_width = unit.c(unit(2, "mm"), max_text_width(txt)))
    } else right_anno <- NULL
    
    # 2nd heatmap:
    hm2 <- switch(type, 
        DA = {
            # relative cluster abundances by sample
            ns <- table(x$cluster_id, x$sample_id)
            fq <- prop.table(ns, 1)
            fq <- fq[top$cluster_id, ]
            fq <- as.matrix(unclass(fq))
            if (normalize)
                fq <- .z_normalize(asin(sqrt(fq)))
            Heatmap(
                matrix = fq,
                name = paste0("normalized\n"[normalize], "frequency"),
                col = hm2_pal,
                row_title = "cluster_id",
                column_title = "sample_id",
                cluster_rows = !order, 
                cluster_columns = FALSE,
                show_row_names = is.null(hm1), 
                row_names_side = "left",
                column_title_side = "bottom",
                top_annotation = col_anno,
                #row_names_gp = gpar(fontsize = 8),
                #column_names_gp = gpar(fontsize = 8),
                rect_gp = gpar(col = "white"),
                right_annotation = right_anno)
        },
        DS = {
            # median state-marker expression by sample
            es <- assay(x, "exprs")
            cs_by_ks <- .split_cells(x, c("cluster_id", "sample_id"))
            ms_by_ks <- t(mapply(function(k, g)
                vapply(cs_by_ks[[k]], function(cs) {
                    if (length(cs) == 0)
                        return(NA)
                    get(fun)(es[g, cs, drop=FALSE])
                }, numeric(1)),
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
            Heatmap(
                matrix = ms_by_ks,
                name = paste0("z-normalized\n"[normalize], "expression"),
                col = hm2_pal,
                column_title = "sample_id",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                top_annotation = col_anno,
                column_title_side = "bottom",
                show_row_names = is.null(hm1) || is.null(right_anno),
                row_names_side = ifelse(is.null(hm1), "left", "right"),
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "median",
                #row_names_gp = gpar(fontsize = 8),
                #column_names_gp = gpar(fontsize = 8),
                rect_gp = gpar(col = "white"),
                right_annotation = right_anno)
        })
    hm1 + hm2
}
