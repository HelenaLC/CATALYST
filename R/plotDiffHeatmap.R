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
#' @param k character string specifying 
#'   the clustering in \code{x} from which \code{y} was obtained.
#'   If NULL, \code{plotDiffHeatmap} will try and guess it,
#'   which will be inaccurate if multiple clusterings share the same levels.
#' @param top_n numeric. Number of top clusters (if \code{type = "DA"})
#'   or cluster-marker combinations (if \code{type = "DS"}) to display.
#' @param fdr numeric threshold on adjusted p-values below which 
#'   results should be retained and considered to be significant.
#' @param lfc numeric threshold on logFCs above which to retain results.
#' @param all logical specifying whether all \code{top_n} results should 
#'   be displayed. If \code{TRUE}, \code{fdr,lfc} filtering is skipped.
#' @param sort_by character string specifying the \code{y} column to sort by; 
#'   \code{"none"} to retain original ordering. Adj. p-values will increase, 
#'   logFCs will decreasing from top to bottom.
#' @param y_cols named list specifying columns in \code{y} that contain
#'   adjusted p-values (\code{padj}), logFCs (\code{lfc}) and, 
#'   for DS results, feature names (\code{target}).
#'   When only some \code{y_cols} differ from the defaults,
#'   specifying only these is sufficient.
#' @param assay character string specifying which assay 
#'   data to use; valid values are \code{assayNames(x)}.
#' @param fun character string specifying the function to use 
#'   as summary statistic for aggregation of \code{assay} data.
#' @param normalize logical specifying whether Z-score normalized values 
#'   should be plotted. If \code{y} contains DA analysis results, 
#'   frequencies will be arcsine-square-root scaled prior to normalization.
#' @param row_anno logical specifying whether to include a row annotation 
#'   indicating whether cluster (DA) or cluster-marker combinations (DS) 
#'   are significant, labeled with adjusted p-values, as well as logFCs.
#' @param col_anno logical specifying whether to include column annotations 
#'   for all non-numeric cell metadata variables; or a character vector 
#'   in \code{names(colData(x))} to include only a subset of annotations.
#'   (Only variables that map uniquely to each sample will be included)
#' @param hm_pal character vector of colors 
#'   to interpolate for the heatmap. Defaults to \code{brewer.pal}'s 
#'   \code{"RdYlBu"} for DS, \code{"RdBu"} for DA results heatmaps.
#' @param fdr_pal,lfc_pal character vector of colors to use for row annotations
#' \itemize{
#' \item{\code{fdr_pal}}{length 2 for (non-)significant at given \code{fdr}}
#' \item{\code{lfc_pal}}{length 3 for negative, zero and positive}}
#' 
#' @return a \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#' 
#' @author Lukas M Weber & Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce, verbose = FALSE)
#' 
#' ## differential analysis
#' library(diffcyt)
#' 
#' # create design & constrast matrix
#' design <- createDesignMatrix(ei(sce), cols_design=2:3)
#' contrast <- createContrast(c(0, 1, 0, 0, 0))
#' 
#' # test for
#' # - differential abundance (DA) of clusters
#' # - differential states (DS) within clusters
#' 
#' da <- diffcyt(sce, design = design, contrast = contrast, 
#'     analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
#'     clustering_to_use = "meta20", verbose = FALSE)
#'     
#' ds <- diffcyt(sce, design = design, contrast = contrast, 
#'     analysis_type = "DS", method_DS = "diffcyt-DS-limma",
#'     clustering_to_use = "meta20", verbose = FALSE)
#'     
#' # extract result tables
#' da <- rowData(da$res)
#' ds <- rowData(ds$res)
#'     
#' # display test results for
#' # - top DA clusters
#' # - top DS cluster-marker combinations
#' plotDiffHeatmap(sce, da)
#' plotDiffHeatmap(sce, ds)
#' 
#' # visualize results for subset of clusters
#' sub <- filterSCE(sce, cluster_id %in% seq_len(5), k = "meta20")
#' plotDiffHeatmap(sub, da, all = TRUE, sort_by = "none")
#' 
#' # visualize results for selected feature
#' # & include only selected annotation
#' plotDiffHeatmap(sce["pp38", ], ds, col_anno = "condition", all = TRUE)
#' 
#' @importFrom ComplexHeatmap rowAnnotation row_anno_text Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom dplyr rename mutate_if
#' @importFrom grid gpar unit
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales scientific
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom S4Vectors metadata
#' @export

plotDiffHeatmap <- function(x, y, k = NULL,
    top_n = 20, fdr = 0.05, lfc = 1, all = FALSE,
    sort_by = c("padj", "lfc", "none"), 
    y_cols = list(padj = "p_adj", lfc = "logFC", target = "marker_id"),
    assay = "exprs", fun = c("median", "mean", "sum"), 
    normalize = TRUE, col_anno = TRUE, row_anno = TRUE,
    hm_pal = NULL, 
    fdr_pal = c("lightgrey", "lightgreen"),
    lfc_pal = c("blue3", "white", "red3")) {
    
    # check validity of input arguments
    fun <- match.arg(fun)
    sort_by <- match.arg(sort_by)
    args <- as.list(environment())
    
    defs <- as.list(formals("plotDiffHeatmap")$y_cols[-1])
    miss <- !names(defs) %in% names(args$y_cols)
    if (any(miss)) y_cols <- args$y_cols <- 
        c(args$y_cols, defs[miss])[names(defs)]
    
    .check_args_plotDiffHeatmap(args)
    stopifnot(y_cols[[sort_by]] %in% names(y))
    y_cols <- y_cols[y_cols %in% names(y)]
    
    # guess clustering to use
    if (is.null(k)) {
        kids <- levels(y$cluster_id)
        same <- vapply(cluster_codes(x), function(u) 
            identical(levels(u), kids), logical(1))
        if (!any(same)) 
            stop("Couldn't match any clustering",
                " in input data 'x' with results in 'y'.")
        k <- names(cluster_codes(x))[same][1]
    } else {
        k <- .check_k(x, k)
    }
    x$cluster_id <- cluster_ids(x, k)
    
    # get feature column
    y <- data.frame(y, check.names = FALSE)
    y <- mutate_if(y, is.factor, as.character)
    if (any(rownames(x) %in% unlist(y))) {
        features <- intersect(rownames(x), y[[y_cols$target]])
        if (length(features) == 0)
            stop("Couldn't match features between",
                " results 'y' and input data 'x'.")
        i <- y[[y_cols$target]] %in% rownames(x)
        type <- "ds"
    } else {
        i <- TRUE
        type <- "da"
    }

    # rename relevant result variables
    y <- rename(y, 
        target = y_cols$target,
        padj = y_cols$padj, 
        lfc = y_cols$lfc)
    
    # filter results
    i <- i & !is.na(y$padj) & y$cluster_id %in% levels(x$cluster_id)
    if (!all) {
        i <- i & y$padj < fdr
        if (!is.null(y$lfc))
            i <- i & abs(y$lfc) > lfc
    }
    y <- y[i, , drop = FALSE]

    if (nrow(y) == 0)
        stop("No results remaining;",
            " perhaps 'x' or 'y' has been filtered,",
            " or features couldn't be matched.")
    
    # get clusters/cluster-marker combinations to plot
    if (sort_by != "none") {
        o <- order(abs(y[[sort_by]]), 
            decreasing = (sort_by == "lfc"))
        y <- y[o, , drop = FALSE]
    }
    if (top_n > nrow(y)) 
        top_n <- nrow(y)
    top <- y[seq_len(top_n), ]
    
    # column annotation of non-numeric cell metadata variables
    if (!isFALSE(col_anno)) {
        sids <- levels(droplevels(factor(x$sample_id)))
        top_anno <- .anno_factors(x, sids, col_anno, "column")
    } else top_anno <- NULL
    
    if (is.null(hm_pal)) hm_pal <- rev(brewer.pal(11, 
        ifelse(type == "ds", "RdYlBu", "RdBu")))
    
    # row annotation: significant = (adj. p-values <= th)
    if (row_anno) {
        s <- factor(
            ifelse(top$padj < fdr, "yes", "no"), 
            levels = c("no", "yes"))
        if (!is.null(top$lfc)) {
            lfc_lims <- range(top$lfc, na.rm = TRUE)
            if (all(lfc_lims > 0)) {
                lfc_brks <- c(0, lfc_lims[2])
                lfc_pal <- lfc_pal[-1]
            } else if (all(lfc_lims < 0)) {
                lfc_brks <- c(lfc_lims[1], 0)
                lfc_pal <- lfc_pal[-3]
            } else lfc_brks <- c(lfc_lims[1], 0, lfc_lims[2])
            lfc_anno <- top$lfc
            anno_cols <- list(logFC = colorRamp2(lfc_brks, lfc_pal))
        } else {
            lfc_anno <- NULL
            anno_cols <- list()
        }
        names(fdr_pal) <- levels(s)
        anno_cols$significant <- fdr_pal
        right_anno <- rowAnnotation(
            logFC = lfc_anno,
            significant = s,
            "foo" = row_anno_text(
                scientific(top$padj, 2),
                gp = gpar(fontsize = 8)),
            col = anno_cols,
            gp = gpar(col = "white"),
            show_annotation_name = FALSE,
            simple_anno_size = unit(4, "mm"))
    } else right_anno <- NULL

    switch(type, 
        # relative cluster abundances by sample
        da = {
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
        ds = {
            y <- assay(x, assay)
            cs <- .split_cells(x, c("cluster_id", "sample_id"))
            z <- t(mapply(function(k, g)
                vapply(cs[[k]], function(cs) {
                    if (length(cs) == 0) return(NA)
                    get(fun)(y[g, cs, drop = FALSE])
                }, numeric(1)),
                k = top$cluster_id, 
                g = top$target))
            rownames(z) <- sprintf("%s(%s)", top$target, top$cluster_id)
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
                right_annotation = right_anno,
                heatmap_legend_param = list(title_gp = gpar(
                    fontsize = 10, fontface = "bold", lineheight = 0.8)))
    }) 
}
