#' @rdname clrDR
#' @title DR plot on CLR of proportions
#'
#' @description
#' Computes centered log-ratios (CLR) on cluster/sample proportions across
#' samples/clusters, and visualizes them in a lower-dimensional space,
#' highlighting differences in composition between samples/clusters.
#'
#' @param by character string specifying across which IDs to compute CLRs
#' \itemize{
#' \item{\code{by = "sample_id"} compute CLRs across
#'   relative abundances of samples across clusters;
#'   each point in the embedded space represents a sample.}
#' \item{\code{by = "cluster_id"} compute CLRs across
#'   relative abundances of clusters across samples;
#'   each point in the embedded space represents a cluster.}
#' }
#' @param k character string specifying which clustering to use;
#'   valid values are \code{names(cluster_codes(x))}.
#' @param base integer scalar specifying the logarithm base to use.
#' @param arrows logical specifying whether to include arrows for PC loadings.
#' @param point_col,arrow_col character string specifying a non-numeric
#'   cell metadata column to color points and PC loading arrows by;
#'   valid values are \code{names(colData(x))}.
#' @param arrow_len non-zero single numeric specifying the length of loading
#'   vectors relative to the largest xy-coordinate in the embedded space;
#'   NULL for no re-sizing (see details).
#' @param arrow_opa single numeric in [0,1] specifying the opacity (alpha)
#'   of PC loading arrows when they are grouped; 0 will hide individual arrows.
#' @param label_by character string specifying a non-numeric sample metadata
#'   variable to label points by; valid values are \code{names(colData(x))}.
#' @param size_by logical specifying whether to scale point sizes by the number
#'   of cells in a given sample/cluster (for \code{by = "sample/cluster_id"}).
#' @param point_pal,arrow_pal character string of colors to use
#'   for points and PC loading arrows. Arguments default to
#'   \code{.cluster_cols} for clusters (defined internally),
#'   and \code{brewer.pal}'s \code{"Set3"} for samples.
#' @inheritParams runDR
#' @inheritParams pbMDS
#' @inheritParams plotDR
#'
#' @details
#' \describe{
#' \item{The centered log-ratio (CLR)}{
#' Let \code{k} be one of \eqn{S} samples, \code{k} one of \eqn{K} clusters,
#' and \code{p(s,k)} be the proportion of cells from \code{s} in \eqn{k}.
#' The centered log-ratio (CLR) is defined as
#' \deqn{clr(sk) = log p(s,k) - \sum p(s,k) / K}
#' and analogous for clusters replacing \code{s} by \code{k} and \code{K} by
#' \code{S}. Thus, each sample/cluster gives a vector with length \code{K/S}
#' and mean 0, and the CLRs computed across all instances can be represented
#' as a matrix with dimensions \code{S} x \code{K} (or \code{K} x \code{S}
#' for clusters) that we embed into a lower dimensional space.}
#'
#' \item{Dimensionality reduction}{
#' In principle, \code{clrDR} allows any dimension reduction to be applied on
#' the CLRs. The default method (\code{dr = "PCA"}) will include the percentage
#' of variance explained by each principal component (PC) in the axis labels.
#'
#' Noteworthily, distances between points in the lower-dimensional space are
#' meaningful only for linear DR methods (PCA and MDS), and results obtained
#' from other methods should be interpreted with caution. Thus, the output
#' plot's aspect ratio should be kept as is for PCA and MDS; non-linear
#' DR methods can use \code{aspect.ratio = 1}, rendering a square plot.}
#'
#' \item{Interpreting PC loadings}{
#' For \code{dr = "PCA"}, PC loadings will be represented as arrows that may be
#' interpreted as follows: 0° (180°) between vectors indicates a strong positive
#' (negative) relation between them, while vectors that are orthogonal to each
#' another (90°) are roughly independent.
#'
#' When a vector points towards a given quadrant, the variability in proportions
#' for the points within this quadrant are largely driven by the corresponding
#' variable. Here, only the relative orientation of vectors to one another and
#' to the PC axes is meaningful; however, the sign of loadings (i.e., whether
#' an arrow points left or right) can be flipped when re-computing PCs.
#'
#' When \code{arrow_len} is specified, PC loading vectors will be re-scaled to
#' improve their visibility. Here, a value of 1 will stretch vectors such that
#' the largest loading will touch on the outer most point. Importantly, while
#' absolute arrow lengths are not interpretable, their relative length is.}}
#'
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#'
#' @return a \code{ggplot} object.
#'
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce)
#'
#' # CLR on sample proportions across clusters
#' # (1st vs. 3rd PCA; include sample labels)
#' clrDR(sce, by = "sample_id", k = "meta12",
#'   dims = c(1, 3), label_by = "sample_id")
#'
#' # CLR on cluster proportions across samples
#' # (use custom colors for both points & loadings)
#' clrDR(sce, by = "cluster_id",
#'   point_pal = hcl.colors(10, "Spectral"),
#'   arrow_pal = c("royalblue", "orange"))
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SummarizedExperiment colData
#' @importFrom scater calculatePCA calculateMDS
#'   calculateUMAP calculateTSNE calculateDiffusionMap
#' @export

clrDR <- function(x,
    dr = c("PCA", "MDS", "UMAP", "TSNE", "DiffusionMap"),
    by = c("sample_id", "cluster_id"), k = "meta20",
    dims = c(1, 2), base = 2, arrows = TRUE,
    point_col = switch(by, sample_id = "condition", "cluster_id"),
    arrow_col = switch(by, sample_id = "cluster_id", "condition"),
    arrow_len = 0.5, arrow_opa = 0.5,
    label_by = NULL, size_by = TRUE,
    point_pal = NULL, arrow_pal = NULL) {

    # check validity of input arguments
    by <- match.arg(by)
    dr <- match.arg(dr)
    args <- as.list(environment())
    .check_args_clrDR(args)

    # if unspecified, get default color palettes
    k_pal <- .cluster_cols
    s_pal <- brewer.pal(8, "Set3")[-2]
    if (is.null(point_pal)) point_pal <- switch(by, cluster_id = k_pal, s_pal)
    if (is.null(arrow_pal)) arrow_pal <- switch(by, sample_id = k_pal, s_pal)

    # compute CLR on cluster-proportions across samples
    x$cluster_id <- cluster_ids(x, k)
    ncs <- table(x$sample_id, x$cluster_id)
    fqs <- log(prop.table(ncs + 1, 1), base)
    clr <- fqs-rowMeans(fqs)
    if (by == "sample_id") clr <- t(clr)
    #clr <- scale(t(clr))

    # run dimensionality reduction
    fun <- get(paste0("calculate", dr))
    xy <- fun(clr, ncomponents = max(dims))

    # construct data.frame for plotting
    df <- data.frame(xy[, dims])
    colnames(df) <- c("x", "y")
    df[[by]] <- factor(colnames(clr), levels(x[[by]]))

    # include available sample metadata
    if (!isTRUE(by == "cluster_id")) {
        m <- match(df$sample_id, x$sample_id)
        i <- setdiff(names(colData(x)), names(df))
        df <- cbind(df, colData(x)[m, i])
    }

    # add sample cell counts
    if (size_by) {
        size_by <- "n_cells"
        df$n_cells <- tabulate(x[[by]])
    } else size_by <- NULL

    # get axis labels & number of legend columns to use
    ncol <- ifelse(!is.null(point_col) && nlevels(df[[point_col]]) > 10, 2, 1)
    labs <- paste(switch(dr, PCA = "PC", paste(dr, "dim. ")), dims)

    # ramp point color palette
    np <- nlevels(df[[point_col]])
    if (length(point_pal) < np)
        point_pal <- colorRampPalette(point_pal)(np)
    pals <- list(p = point_pal, a = arrow_pal)

    if (dr %in% c("PCA", "MDS")) {
        lines <- list(
            geom_vline(xintercept = 0, lty = 2),
            geom_hline(yintercept = 0, lty = 2))
        asp <- coord_equal()
    } else lines <- asp <- NULL

    # make base plot
    p <- ggplot(df, aes_string("x", "y", fill = point_col)) +
        lines + guides(
            fill = guide_legend(order = 1, ncol = ncol,
                override.aes = list(alpha = 1, size = 4)),
            size = guide_legend(order = 2,
                override.aes = list(shape = 19))) +
        (if (!is.null(point_pal))
            scale_fill_manual(
                values = point_pal)) +
        (if (!is.null(label_by))
            geom_text_repel(
                show.legend = FALSE,
                aes_string(label = label_by))) +
        (if (is.null(size_by)) {
            geom_point(size = 5, alpha = 0.8, shape = 21, stroke = 0)
        } else geom_point(
            aes_string(size = size_by),
            alpha = 0.8, shape = 21, stroke = 0)) +
        labs(x = labs[1], y = labs[2]) +
        asp + theme_minimal() + theme(
            panel.grid.minor = element_blank(),
            legend.key.height  =  unit(0.8, "lines"),
            axis.text = element_text(color = "black"),
            aspect.ratio = if (is.null(asp)) 1 else NULL)

    if (dr == "PCA") {
        # include variance explained in axis labels
        ve <- attr(xy, "percentVar")[dims]
        ve <- sprintf("(%s%%)", round(ve))
        labs <- unlist(p$labels[c("x", "y")])
        p$labels[c("x", "y")] <- paste(labs, ve)
        # (optionally) add loading arrows
        if (arrows) {
            # construct data.frame of PC loadings
            rot <- attr(xy, "rotation")
            rot <- data.frame(rot)
            colnames(rot) <- c("x", "y")
            # add sample metadata
            arrow_by <- switch(by,
                sample_id = "cluster_id",
                cluster_id = "sample_id")
            if (arrow_by == "sample_id") {
                m <- match(rownames(rot), x$sample_id)
                rot <- cbind(rot, colData(x)[m, ])
            }
            rot[[arrow_by]] <- factor(rownames(rot), levels(x[[arrow_by]]))
            # scale arrow lengths
            if (!is.null(arrow_len)) {
                xy <- rot[c("x", "y")]
                m1 <- max(abs(df[c("x", "y")]))
                m2 <- max(abs(xy))
                rot[c("x", "y")] <- m1/m2*arrow_len*xy
            }
            # ramp arrow color palette
            na <- nlevels(x[[arrow_col]])
            if (length(arrow_pal) < na)
                arrow_pal <- colorRampPalette(arrow_pal)(na)
            p <- p + geom_segment(
                data = rot, inherit.aes = FALSE,
                linewidth = 0.5, arrow = arrow(length = unit(1, "mm")),
                aes_string(0, 0, xend = "x", yend = "y", col = arrow_col)) +
                scale_color_manual(values = arrow_pal) +
                guides(col = guide_legend(
                    override.aes = list(size = 1),
                    ncol = ifelse(nrow(rot) > 10, 2, 1)))
        }
    }
    return(p)
}
