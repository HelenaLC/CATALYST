#' @rdname plotExprHeatmap
#' @title Plot expression heatmap
#' 
#' @description 
#' Heatmap of marker expressions aggregated by sample, cluster, or both;
#' with options to include annotation of cell metadata factors, clustering(s),
#' as well as relative and absolute cell counts.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param features character string specifying which features to include;
#'   valid values are \code{"type"/"state"} for \code{type/state_markers(x)} 
#'   if \code{rowData(x)$marker_class} have been specified; 
#'   a subset of \code{rownames(x)}; NULL to use all features.
#'   When \code{by = "both"}, only 1 feature is allowed.
#' @param by character string specifying 
#'   whether to aggregate by sample, cluster, both.
#' @param k character string specifying which 
#'   clustering to use when \code{by != "sample_id"};
#'   \code{assay} data will be aggregated across these cluster IDs.
#' @param m character string specifying a metaclustering to include as an 
#'   annotation when \code{by != "sample_id"} and \code{row_anno = TRUE}.
#' @param assay character string specifying which assay 
#'   data to use; valid values are \code{assayNames(x)}.
#' @param fun character string specifying 
#'   the function to use as summary statistic.
#' @param scale character string specifying the scaling strategy:
#' \itemize{
#'   \item{\code{"first"}: scale & trim then aggregate}
#'   \item{\code{"last"}: aggregate then scale & trim}
#'   \item{\code{"never"}: aggregate only}
#' } If \code{scale != "never"}, data will be scaled using lower 
#'   (\code{q}\%) and upper (\code{1-q}\%) quantiles as boundaries.
#' @param q single numeric in [0,0.5) determining the 
#'   quantiles to trim when \code{scale != "never"}.
#' @param row_anno,col_anno logical specifying whether to include row/column 
#'   annotations (see details); when one axis corresponds to samples 
#'   (\code{by != "cluster_id"}), this can be a character vector specifying 
#'   a subset of \code{names(colData(x))} to be included as annotations.
#' @param row_clust,col_clust logical specifying whether rows/columns 
#'   should be hierarchically clustered and re-ordered accordingly.
#' @param row_dend,col_dend logical specifying 
#'   whether to include the row/column dendrograms.
#' @param bars logical specifying whether to include a barplot 
#'   of cell counts per cluster as a right-hand side row annotation.
#' @param perc logical specifying whether to display 
#'   percentage labels next to bars when \code{bars = TRUE}.
#' @param bin_anno logical specifying whether to display values inside bins.
#' @param hm_pal character vector of colors to interpolate for the heatmap. 
#' @param k_pal,m_pal character vector of colors to interpolate 
#'   for cluster annotations when \code{by != "sample_id"}.
#' @param distance character string specifying the distance metric 
#'  to use for both row and column hierarchical clustering; 
#'  passed to \code{\link[ComplexHeatmap]{Heatmap}} 
#' @param linkage character string specifying the agglomeration method 
#'  to use for both row and column hierarchical clustering; 
#'  passed to \code{\link[ComplexHeatmap]{Heatmap}} 
#' 
#' @return a \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#' 
#' @details 
#'   By default (\code{row/col_anno = TRUE}), for axes corresponding to samples 
#'   (y-axis for \code{by = "sample_id"} and x-axis for \code{by = "both"}), 
#'   annotations will be drawn for all non-numeric cell metadata variables.
#'   Alternatively, a specific subset of annotations can be included
#'   for only a subset of variables by specifying \code{row/col_anno} 
#'   to be a character vector in \code{names(colData(x))} (see examples).  
#'   
#'   For axes corresponding to clusters (y-axis for \code{by = "cluster_id"} 
#'   and \code{"both"}), annotations will be drawn for the specified 
#'   clustering(s) (arguments \code{k} and \code{m}).
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
#' \code{\link{plotFreqHeatmap}}, 
#' \code{\link{plotMultiHeatmap}}
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce) 
#' 
#' # median scaled & trimmed expression by cluster
#' plotExprHeatmap(sce, 
#'   by = "cluster_id", k = "meta8",
#'   scale = "first", q = 0.05, bars = FALSE)
#' 
#' # scale each marker between 0 and 1 
#' # after aggregation (without trimming)
#' plotExprHeatmap(sce, 
#'   scale = "last", q = 0,
#'   bars = TRUE, perc = TRUE,
#'   hm_pal = hcl.colors(10, "YlGnBu", rev = TRUE))
#' 
#' # raw (un-scaled) median expression by cluster-sample
#' plotExprHeatmap(sce,
#'   features = "pp38", by = "both", k = "meta10", 
#'   scale = "never", row_anno = FALSE, bars = FALSE)
#'   
#' # include only subset of samples
#' sub <- filterSCE(sce, 
#'   patient_id != "Patient",
#'   sample_id != "Ref3")
#'  
#' # includes specific annotations &
#' # split into CDx & all other markers
#' is_cd <- grepl("CD", rownames(sce))
#' plotExprHeatmap(sub, 
#'   rownames(sce)[is_cd], 
#'   row_anno = "condition", 
#'   bars = FALSE)
#' plotExprHeatmap(sub, 
#'   rownames(sce)[!is_cd], 
#'   row_anno = "patient_id",
#'   bars = FALSE)
#' 
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom dplyr select select_if
#' @importFrom grid gpar grid.text
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay
#' @export

plotExprHeatmap <- function(x, features = NULL, 
    by = c("sample_id", "cluster_id", "both"), k = "meta20", m = NULL,
    assay = "exprs", fun = c("median", "mean", "sum"), 
    scale = c("first", "last", "never"), q = 0.01, 
    row_anno = TRUE, col_anno = TRUE,
    row_clust = TRUE, col_clust = TRUE, 
    row_dend = TRUE, col_dend = TRUE, 
    bars = FALSE, perc = FALSE, bin_anno = FALSE,
    hm_pal = rev(brewer.pal(11, "RdYlBu")), 
    k_pal = CATALYST:::.cluster_cols, m_pal = k_pal,
    distance = c(
        "euclidean", "maximum", "manhattan", 
        "canberra", "binary", "minkowski"), 
    linkage = c(
        "average", "ward.D", "single", "complete", 
        "mcquitty", "median", "centroid", "ward.D2")) {
    
    # check validity of input arguments
    args <- as.list(environment())
    .check_args_plotExprHeatmap(args)
    distance <- match.arg(distance)
    linkage <- match.arg(linkage)
    scale <- match.arg(scale)
    fun <- match.arg(fun)
    by <- match.arg(by)

    # subset features of interest
    x <- x[unique(.get_features(x, features)), ]
    
    # get specified cluster IDs
    if (by != "sample_id") {
        .check_k(x, k)
        x$cluster_id <- cluster_ids(x, k)
    } 
    if (by == "both")
        by <- c("cluster_id", "sample_id")

    # aggregate to pseudobulks by sample/cluster/both
    # using 'assay' data & 'fun' as summary statistic
    .do_agg <- function() {
        z <- .agg(x, by, fun, assay)
        if (length(by) > 1) {
            z <- do.call("rbind", z)
            rownames(z) <- levels(x$cluster_id)
        }
        return(z)
    }
    # do 0-1 scaling for each marker trimming 
    # lower ('q'%) & upper (1-'q'%) quantiles
    .do_scale <- function() {
        if (scale == "first") {
            z <- assay(x, assay)
            z <- .scale_exprs(z, 1, q)
            assay(x, assay, FALSE) <- z
            return(x)
        } else .scale_exprs(z, 1, q)
    }
    
    # apply one of...
    # - scale & trim then aggregate
    # - aggregate then scale & trim
    # - aggregate only
    z <- switch(scale,
        first = { x <- .do_scale(); .do_agg() },
        last =  { z <- .do_agg(); .do_scale() },
        never = { .do_agg() })
    if (length(by) == 1) z <- t(z)
    
    if (scale != "never" && !(assay == "counts" && fun == "sum")) {
        qs <- round(quantile(z, c(0.01, 0.99))*5)/5
        lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
    } else lgd_aes <- list()
    lgd_aes$title_gp <- gpar(
        fontsize = 10, 
        fontface = "bold", 
        lineheight = 0.8)
    
    # left-hand side heatmap annotation:
    # non-numeric cell metadata variables
    sids <- levels(droplevels(factor(x$sample_id)))
    if (!isFALSE(row_anno)) {
        left_anno <- switch(by[1],
            sample_id = .anno_factors(x, sids, row_anno, "row"),
            .anno_clusters(x, k, m, k_pal, m_pal))
    } else left_anno <- NULL
    if (!isFALSE(col_anno) && length(by) == 2) {
        top_anno <- .anno_factors(x, sids, col_anno, "colum")
    } else top_anno <- NULL

    # right-hand side heatmap annotation:
    # labeled barplot of event counts by sample
    if (bars) {
        right_anno <- .anno_counts(x[[by[1]]], perc)
    } else right_anno <- NULL
    
    # get bin annotation
    if (bin_anno) {
        cell_fun <- function(j, i, x, y, ...) 
            grid.text(
                gp = gpar(fontsize = 8), 
                sprintf("%.2f", z[i, j]), x, y)
    } else cell_fun <- NULL

    a <- ifelse(assay == "exprs", "expression", assay)
    f <- switch(fun, "median" = "med", fun)
    hm_title <- switch(scale, 
        first = sprintf("%s %s\n%s", fun, "scaled", a),
        last = sprintf("%s %s\n%s", "scaled", fun, a),
        never = paste(fun, a, sep = "\n"))
    if (length(by) == 2) {
        col_title <- features
    } else if (length(features) == 1 
        && features %in% c("type", "state")) {
        col_title <- paste0(features, "_markers")
    } else col_title <- ""
    
    Heatmap(
        matrix = z,
        name = hm_title,
        col = colorRamp2(
            seq(min(z), max(z), l = n <- 100),
            colorRampPalette(hm_pal)(n)),
        column_title = col_title,
        column_title_side = ifelse(length(by) == 2, "top", "bottom"),
        cell_fun = cell_fun, 
        cluster_rows = row_clust,
        cluster_columns = col_clust,
        show_row_dend = row_dend,
        show_column_dend = col_dend,
        clustering_distance_rows = distance,
        clustering_method_rows = linkage,
        clustering_distance_columns = distance,
        clustering_method_columns = linkage,
        show_row_names = (
            is.null(left_anno) 
            || isTRUE(by == "sample_id")) && !perc,
        row_names_side = ifelse(
            by[1] == "cluster_id"
            || isFALSE(row_anno) && !row_dend 
            || isFALSE(row_clust), 
            "left", "right"),
        top_annotation = top_anno,
        left_annotation = left_anno,
        right_annotation = right_anno,
        rect_gp = gpar(col = "white"),
        heatmap_legend_param = lgd_aes)
}
