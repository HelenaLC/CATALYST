#' @rdname plotClusterHeatmap
#' @title Plot cluster heatmap
#' 
#' @description 
#' Plots heatmaps summarizing a clustering and/or metaclustering of interest.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param hm2 character string. Specifies the right-hand side heatmap. 
#'   One of: \itemize{
#'   \item{\code{"abundances"}: cluster frequencies across samples}
#'   \item{\code{"state"}: median state-marker expressions 
#'     across clusters (analogous to the left-hand side heatmap)}
#'   \item{a character string/vector corresponding to one/multiple marker(s): 
#'     median marker expressions across samples and clusters}}
#' @param k character string specifying the clustering 
#'   across which median marker expressions should be computed;
#'   valid values are \code{names(cluster_codes(x))}.
#' @param m character string specifying a metaclustering to include 
#'   as a row annotation (for display only, does not effect computations!);
#'   valid values are \code{names(cluster_codes(x))}.
#' @param assay character string specifying which assay data to use;
#'   valid values are \code{assayNames(x)}.
#' @param fun character string specifying 
#'   the function to use as summary statistic.
#' @param row_anno logical. Should row annotations for each
#'   clustering specified by \code{k} (and \code{m}) be included?
#' @param split_by character string specifying a cell metadata variable
#'   to group the data by; valid values are \code{names(colData(x))}. 
#'   If specified, a list of multiple heatmaps will be rendered.
#' @param scale logical specifying whether expression values should be scaled
#'   between 0 and 1 using lower (1\%) and upper (99\%) quantiles as boundaries.
#' @param row_dend logical. Should a row dendrogram
#'   for the hierarchical clustering of clusters be included?
#' @param col_dend logical. Should a column dendrogram
#'   for the hierarchical clustering of markers be included?
#' @param draw_freqs logical specifying whether to display
#'   a barplot of cell counts labeled with proportions 
#'   for each cluster as a right-hand side row annotation.
#' @param hm1_pal,hm2_pal character vector of colors 
#'   to interpolate for each heatmap.  
#' @param k_pal,m_pal character vector of colors
#'   to use for cluster and merging row annotations.
#'   If less than \code{nlevels(cluster_ids(x, k/m))} 
#'   values are supplied, colors will be interpolated 
#'   via \code{\link[grDevices]{colorRampPalette}}.
#' 
#' @details 
#' In its 1st panel, \code{plotClusterHeatmap} will display (scaled) 
#' type-marker expressions aggregated by cluster (across all samples).
#' Depending on argument \code{hm2}, the 2nd panel will contain one of:
#' \describe{
#' \item{\code{hm2 = "abundances"}}{
#'   relataive cluster abundances by cluster & sample}
#' \item{\code{hm2 = "state"}}{
#'   aggregated (scaled) state-marker expressions by 
#'   cluster (across all samples; analogous to panel 1)}
#' \item{\code{hm2 \%in\% rownames(x)}}{
#'   aggregated (scaled) marker expression by cluster & sample}
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
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#' 
#' plotClusterHeatmap(sce, hm2="abundances", draw_freqs = TRUE)
#' plotClusterHeatmap(sce, hm2="pS6", k="meta12", m="meta8")
#' plotClusterHeatmap(sce, hm2="state", split_by="condition")
#' 
#' @import ComplexHeatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assay
#' @export

plotClusterHeatmap <- function(x, hm2 = NULL, 
    k = "meta20", m = NULL, 
    assay = "exprs", fun = c("median", "mean", "sum"), 
    row_anno = TRUE, split_by = NULL, scale = TRUE, 
    row_dend = TRUE, col_dend = FALSE, draw_freqs = FALSE, 
    hm1_pal = rev(brewer.pal(11, "RdBu")), 
    hm2_pal = (
        if (!is.null(hm2)) 
            if (isTRUE(hm2 == "abundances")) {
                rev(brewer.pal(11, "PuOr")) 
            } else hm1_pal),
    k_pal = CATALYST:::.cluster_cols, m_pal = k_pal) {

    # check validity of input arguments
    args <- as.list(environment())
    .check_args_plotClusterHeatmap(args)
    fun <- match.arg(fun)

    # ramp color palettes
    hm1_pal <- colorRampPalette(hm1_pal)(100)
    if (!is.null(hm2_pal)) 
        hm2_pal <- colorRampPalette(hm2_pal)(100)
    
    # clustering row annotations
    if (row_anno && !is.null(c(k, m))) {
        left_anno <- .get_row_anno(x, k, m, k_pal, m_pal)
    } else left_anno <- NULL
    
    # hierarchical clustering on cell-type marker medians by cluster
    x$cluster_id <- cluster_ids(x, k)
    ms_by_k <- t(.agg(x, "cluster_id", fun, assay))
    d <- dist(ms_by_k[, type_markers(x)])
    row_clustering <- hclust(d, method = "average")
    
    # split cell indices by colData factor
    cs <- seq_len(ncol(x))
    if (!is.null(split_by)) {
        groups <- split(cs, x[[split_by]]) 
    } else groups <- list(cs)   
    
    # optionally scale expression matrix
    if (scale) {
        y <- assay(x, assay)
        y <- .scale_exprs(y, 1)
        assay(x, "scaled", FALSE) <- y
    }
    hms <- lapply(seq_along(groups), function(i) {
        idx <- groups[[i]]
        cs_by_k <- split(idx, x$cluster_id[idx])
        # left-hand side heatmap -----------------------------------------------
        if (!is.null(split_by)) {
            # aggregate subsetted data
            hm1_es <- t(.agg(x[, idx], "cluster_id", fun, assay)) 
        } else {
            if (scale) {
                # re-aggregate scaled data
                hm1_es <- t(.agg(x, "cluster_id", fun, "scaled"))
            } else {
                # use unscaled data from above
                hm1_es <- ms_by_k
            }
        }
        # right-hand side row annotaion of cluster cell counts
        if (draw_freqs) {
            ncs <- tabulate(x$cluster_id[idx])
            frq <- round(ncs/length(idx)*100, 2)
            txt <- sprintf("%s (%s%%)", levels(x$cluster_id), frq)
            right_anno <- rowAnnotation(
                "Frequency [%]" = row_anno_barplot(
                    x = frq, axis = TRUE, border = FALSE,
                    bar_width = 0.8, width = unit(2, "cm"),
                    gp = gpar(fill = "grey50", col = "white")),
                "foo" = row_anno_text(txt, width = max_text_width(txt)))
        } else right_anno <- NULL
        # combine row (cluster) annotations, heatmap of aggregated 
        # type-marker expression by cluster, cell count bars & labels 
        p <- Heatmap(
            matrix = hm1_es[, type_markers(x)], 
            col = hm1_pal, 
            name = paste0("scaled\n"[scale], 
                ifelse(assay == "exprs", "expression", assay)), 
            rect_gp = gpar(col='white'), 
            na_col="lightgrey", 
            cluster_rows = row_clustering, 
            show_row_dend = row_dend, 
            show_column_dend = col_dend,
            column_title = names(groups)[i][!is.null(split_by)],
            left_annotation = left_anno,
            right_annotation = right_anno)
        
        # right-hand side heatmap ----------------------------------------------
        if (!is.null(hm2)) {
            if (isTRUE(hm2 == "abundances")) {
                # relative cluster abundances by samples
                cs <- table(x$cluster_id[idx], x$sample_id[idx])
                fq <- as.matrix(unclass(prop.table(cs, 1)))
                fq <- fq[, !is.na(colSums(fq)), drop = FALSE]
                p <- p + Heatmap(
                    matrix = fq, 
                    name="frequency",
                    col = hm2_pal,
                    na_col="lightgrey", 
                    rect_gp = gpar(col="white"), 
                    show_row_names = FALSE, 
                    cluster_rows = row_clustering, 
                    cluster_columns = FALSE)
            } else if (isTRUE(hm2 == "state")) {
                # aggregated state-marker expression by cluster
                p <- p + Heatmap(
                    matrix = hm1_es[, state_markers(x)], 
                    col = hm2_pal, 
                    na_col="lightgrey", 
                    rect_gp = gpar(col='white'), 
                    show_heatmap_legend = FALSE, 
                    cluster_rows = row_clustering, 
                    cluster_columns = FALSE)
            } else {
                for (ch in hm2) {
                # aggregated marker expression by samples & clusters
                ms <- .agg(x[ch, idx], 
                    by = c("cluster_id", "sample_id"), fun, 
                    assay = ifelse(scale, "scaled", assay))
                ms <- do.call("rbind", ms)
                rownames(ms) <- levels(x$cluster_id)
                p <- p + Heatmap(
                    matrix = ms, 
                    col = hm2_pal, 
                    na_col = "lightgrey", 
                    cluster_rows = row_clustering, 
                    cluster_columns = FALSE,
                    show_row_names = FALSE,
                    column_title = ch, 
                    show_heatmap_legend = FALSE, 
                    rect_gp = gpar(col = "white"))
                }
            }
        }
        return(p)
    })
    if (is.null(split_by)) hms[[1]] else hms
}
