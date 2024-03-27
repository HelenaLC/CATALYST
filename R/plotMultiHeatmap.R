#' @rdname plotMultiHeatmap
#' @title Multi-panel expression & frequency heatmaps
#' 
#' @description Combines expression and frequency heatmaps from 
#' \code{\link{plotExprHeatmap}} and \code{\link{plotFreqHeatmap}}, 
#' respectively, into a \code{\link[ComplexHeatmap]{HeatmapList}}.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param hm1 character string specifying 
#'   which features to include in the 1st heatmap;
#'   valid values are \code{"type"/"state"} for \code{type/state_markers(x)} 
#'   if \code{rowData(x)$marker_class} have been specified; 
#'   a subset of \code{rownames(x)}; NULL to use all features;
#'   and FALSE to omit the 1st heatmap altogether.
#' @param hm2 character string. Specifies the right-hand side heatmap. 
#'   One of: \itemize{
#'   \item{\code{"abundances"}: cluster frequencies across samples}
#'   \item{\code{"state"}: median state-marker expressions 
#'     across clusters (analogous to the left-hand side heatmap)}
#'   \item{a character string/vector corresponding to one/multiple marker(s): 
#'     median marker expressions across samples and clusters}}
#' @param k character string specifying which;
#'   valid values are \code{names(cluster_codes(x))}.
#' @param m character string specifying a metaclustering 
#'   to include as an annotation when \code{row_anno = TRUE}.
#' @param assay character string specifying which assay 
#'   data to use; valid values are \code{assayNames(x)}.
#' @param fun character string specifying 
#'   the function to use as summary statistic.
#' @param scale character string specifying the scaling strategy; 
#'   for expression heatmaps (see \code{\link{plotExprHeatmap}}).
#' @param q single numeric in [0,1) determining the 
#'   quantiles to trim when \code{scale != "never"}.
#' @param normalize logical specifying whether to Z-score normalize 
#'   cluster frequencies across samples; see \code{\link{plotFreqHeatmap}}.
#' @param row_anno,col_anno logical specifying whether to include 
#'   row/column annotations for cell metadata variables and clustering(s); 
#'   see \code{\link{plotExprHeatmap}} and \code{\link{plotFreqHeatmap}}.
#' @param row_clust,col_clust logical specifying whether rows/columns 
#'   should be hierarchically clustered and re-ordered accordingly.
#' @param row_dend,col_dend logical specifying 
#'   whether to include the row/column dendrograms.
#' @param bars logical specifying whether to include a barplot 
#'   of cell counts per cluster as a right-hand side row annotation.
#' @param perc logical specifying whether to display 
#'   percentage labels next to bars when \code{bars = TRUE}.
#' @param hm1_pal,hm2_pal character vector of colors 
#'   to interpolate for each heatmap.  
#' @param k_pal,m_pal character vector of colors
#'   to use for cluster and merging row annotations.
#'   If less than \code{nlevels(cluster_ids(x, k/m))} 
#'   values are supplied, colors will be interpolated 
#'   via \code{\link[grDevices:colorRamp]{colorRampPalette}}.
#' @param distance character string specifying the distance metric 
#'   to use in \code{\link[stats]{dist}} for hierarchical clustering. 
#' @param linkage character string specifying the agglomeration method 
#'   to use in \code{\link[stats]{hclust}} for hierarchical clustering. 
#' 
#' @details 
#' In its 1st panel, \code{plotMultiHeatmap} will display (scaled) 
#' type-marker expressions aggregated by cluster (across all samples).
#' Depending on argument \code{hm2}, the 2nd panel will contain one of:
#' \describe{
#' \item{\code{hm2 = "abundances"}}{
#'   relataive cluster abundances by cluster & sample}
#' \item{\code{hm2 = "state"}}{
#'   aggregated (scaled) state-marker expressions by 
#'   cluster (across all samples; analogous to panel 1)}
#' \item{\code{hm2 \%in\% rownames(x)}}{
#'   aggregated (scaled) marker expressions by cluster & sample}
#' }
#' 
#' @return a \code{\link[ComplexHeatmap]{HeatmapList-class}} object.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Crowell HL, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @seealso 
#' \code{\link{plotMedExprs}}, 
#' \code{\link{plotAbundances}}, 
#' \code{\link{plotExprHeatmap}}, 
#' \code{\link{plotFreqHeatmap}}
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' # state-markers + cluster frequencies
#' plotMultiHeatmap(sce, 
#'   hm1 = "state", hm2 = "abundances", 
#'   bars = TRUE, perc = TRUE)
#' 
#' # type-markers + marker of interest
#' plotMultiHeatmap(sce, hm2 = "pp38", k = "meta12", m = "meta8")
#' 
#' # both, type- & state-markers
#' plotMultiHeatmap(sce, hm2 = "state")
#' 
#' # plot markers of interest side-by-side 
#' # without left-hand side heatmap
#' plotMultiHeatmap(sce, k = "meta10", 
#'   hm1 = NULL, hm2 = c("pS6", "pNFkB", "pBtk"), 
#'   row_anno = FALSE, hm2_pal = c("white", "black"))
#' 
#' @import ComplexHeatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay
#' @export

plotMultiHeatmap <- function(x, 
    hm1 = "type", hm2 = "abundances", 
    k = "meta20", m = NULL,
    assay = "exprs", fun = c("median", "mean", "sum"), 
    scale = c("first", ifelse(hm2 == "state", "first", "last")), 
    q = c(0.01, ifelse(hm2 == "state", 0.01, 0)), 
    normalize = TRUE,
    row_anno = TRUE, 
    col_anno = TRUE,
    row_clust = TRUE, 
    col_clust = c(TRUE, hm2 == "state"),
    row_dend = TRUE, 
    col_dend = c(TRUE, hm2 == "state"),
    bars = FALSE, perc = FALSE,
    hm1_pal = rev(brewer.pal(11, "RdYlBu")), 
    hm2_pal = if (isTRUE(hm2 == "abundances")) 
        rev(brewer.pal(11, "PuOr")) else hm1_pal,
    k_pal = .cluster_cols, m_pal = k_pal,
    distance = c(
        "euclidean", "maximum", "manhattan", 
        "canberra", "binary", "minkowski"), 
    linkage = c(
        "average", "ward.D", "single", "complete", 
        "mcquitty", "median", "centroid", "ward.D2")) {
    
    if (!(isFALSE(hm1) || length(hm1) == 1 && hm1 %in% c("type", "state"))) 
        hm1 <- .get_features(x, hm1)
    stopifnot(is.character(hm2), all(hm2 %in% rownames(x)) 
        || length(hm2) == 1 && hm2 %in% c("abundances", "state"))

    choices <- formals("plotExprHeatmap")$scale
    for (i in seq_along(scale)) 
        match.arg(scale[i], eval(choices))
    
    # recycle single-length arguments
    names(args) <- args <- c("scale", "q", "col_clust", "col_dend")
    vals <- lapply(args, function(u) eval(sym(u)))
    for (i in args[vapply(vals, length, numeric(1)) == 1])
        assign(i, rep(vals[[i]], 2))
    
    # type-marker expression by cluster
    if (!isFALSE(hm1)) {
        a <- plotExprHeatmap(x, hm1, by = "cluster_id",
            k, m, assay, fun, scale[1], q[1], row_anno, col_anno, 
            row_clust, col_clust[1], row_dend, col_dend[1],
            bars, perc, bin_anno = FALSE, hm1_pal, k_pal, m_pal)
    } else a <- NULL
    if (isTRUE(hm2 == "abundances")) {
        # cluster frequencies by sample
        b <- plotFreqHeatmap(x, k, m, 
            normalize, row_anno = FALSE, col_anno, 
            row_clust, col_clust[2], row_dend, col_dend[2], 
            bars = FALSE, perc, hm2_pal, k_pal, m_pal)
        if (!isFALSE(col_anno))
            for (j in seq_along(b@top_annotation@anno_list))
                b@top_annotation@anno_list[[j]]@label <- NULL
        c <- a + b
    } else if (isTRUE(hm2 == "state")) {
        # state-marker expression by cluster
        b <- plotExprHeatmap(x, "state", by = "cluster_id",
            k, m, assay, fun, scale[2], q[2], row_anno = FALSE, col_anno, 
            row_clust, col_clust[2], row_dend, col_dend[2],
            bars = FALSE, perc, bin_anno = FALSE, hm2_pal, k_pal, m_pal)
        # assure both legends show if ranges or colors differ
        if (!isFALSE(hm1)) {
            equal_name <- identical(
                a@matrix_color_mapping@name,
                b@matrix_color_mapping@name)
            equal_breaks <- identical(
                a@matrix_legend_param$at, 
                b@matrix_legend_param$at) 
            equal_colors <- identical(
                a@matrix_color_mapping@colors, 
                b@matrix_color_mapping@colors)
            if (!(equal_name && equal_breaks && equal_colors))
                b@matrix_color_mapping@name <- paste0(b@name, " ")
            b@name <- "foo"
        }
        c <- a + b
    } else {
        # marker expression by cluster-sample
        first <- TRUE
        for (i in hm2) {
            b <- plotExprHeatmap(x, i, by = "both",
                k, m, assay, fun, scale[2], q[2], row_anno = FALSE, col_anno, 
                row_clust, col_clust[2], row_dend, col_dend[2],
                bars = FALSE, perc, bin_anno = FALSE, hm2_pal, k_pal, m_pal)
            # make heatmap identifiers unique &
            # assure both legends show if ranges differ
            b@name <- b@matrix_color_mapping@name <- i
            b@matrix_legend_param$at <- range(
                b@matrix_color_mapping@levels)
            # remove annotation names
            if (!isFALSE(col_anno) && match(i, hm2) < length(hm2)) {
                for (j in seq_along(b@top_annotation@anno_list))
                    b@top_annotation@anno_list[[j]]@name_param$show <- FALSE
            }
            # remove column names if more than 1 & order is consistent
            if (!col_clust[2] && match(i, hm2) > 1)
                b@column_names_param$show <- FALSE
            a <- a + b
        }
        c <- a
    }
    return(c)
}
