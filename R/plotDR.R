#' @rdname plotDR
#' @title Plot reduced dimensions
#' 
#' @description Dimension reduction plot colored 
#' by expression, cluster, sample or group ID.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param dr character string specifying which dimension reduction to use. 
#'   Should be one of \code{reducedDimNames(x)}; default to the 1st available.
#' @param color_by character string specifying the color coding;
#'   valid values are \code{rownames(sce)} and \code{names(colData(x))}.
#' @param facet_by character string specifying a
#'   non-numeric cell metadata column to facet by; 
#'   valid values are \code{names(colData(x))}.
#' @param ncol integer scalar specifying number of facet columns; 
#'   ignored unless coloring by multiple features without facetting
#'   or coloring by a single feature with facetting.
#' @param assay character string specifying which assay data to use
#'   when coloring by marker(s); valid values are \code{assayNames(x)}.
#' @param scale logical specifying whether \code{assay} data should be scaled
#'   between 0 and 1 using lower (1\%) and upper (99\%) expression quantiles;
#'   ignored if \code{!all(color_by \%in\% rownames(x))}.
#' @param q single numeric in [0,0.5) determining the 
#'   quantiles to trim when \code{scale = TRUE}.
#' @param dims length 2 numeric specifying which dimensions to plot.
#' @param k_pal character string specifying the cluster color palette; 
#'   ignored when \code{color_by} is not one of \code{names(cluster_codes(x))}. 
#'   If less than \code{nlevels(cluster_ids(x, k))} are supplied, colors will 
#'   be interpolated via \code{\link[grDevices:colorRamp]{colorRampPalette}}.
#' @param a_pal character string specifying the \code{assay} data palette 
#'   when coloring by feature(s), i.e. \code{all(color_by \%in\% rownames(x))}.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Crowell HL, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' # construct SCE & run clustering
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run clustering & dimension reduction
#' sce <- cluster(sce)
#' sce <- runDR(sce, dr = "UMAP", cells = 100)
#' 
#' # color by single marker, split by sample
#' plotDR(sce, color_by = "CD7", facet_by = "sample_id", ncol = 4)
#' 
#' # color by a set of markers using custom color palette
#' cdx <- grep("CD", rownames(sce), value = TRUE)
#' plotDR(sce, color_by = cdx, ncol = 4, 
#'   a_pal = rev(hcl.colors(10, "Spectral")))
#' 
#' # color by scaled expression for 
#' # set of markers, split by condition
#' plotDR(sce, 
#'   scale = TRUE, 
#'   facet_by = "condition",
#'   color_by = sample(rownames(sce), 4))
#' 
#' # color by 8 metaclusters using custom 
#' # cluster color palette, split by sample
#' p <- plotDR(sce, 
#'   color_by = "meta8", 
#'   facet_by = "sample_id", 
#'   k_pal = c("lightgrey", "cornflowerblue", "navy")) 
#' p$facet$params$ncol <- 4; p
#' 
#' @import ggplot2
#' @importFrom grDevices colorRampPalette hcl.colors
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom SummarizedExperiment assay colData
#' @importFrom stats reformulate
#' @export

plotDR <- function(x, dr = NULL, 
    color_by = "condition", facet_by = NULL, ncol = NULL,
    assay = "exprs", scale = TRUE, q = 0.01, dims = c(1, 2),
    k_pal = .cluster_cols, a_pal = hcl.colors(10, "Viridis")) {
    
    # check validity of input arguments
    stopifnot(
        is(x, "SingleCellExperiment"),
        .check_assay(x, assay),
        length(reducedDims(x)) != 0,
        is.logical(scale), length(scale) == 1,
        is.numeric(q), length(q) == 1, q >= 0, q < 0.5)
    .check_pal(a_pal)
    .check_cd_factor(x, facet_by)
    
    if (!is.null(ncol)) 
        stopifnot(is.numeric(ncol), length(ncol) == 1, ncol %% 1 == 0)
    
    if (is.null(dr)) {
        dr <- reducedDimNames(x)[1]
    } else {
        stopifnot(
            is.character(dr), length(dr) == 1, 
            dr %in% reducedDimNames(x))
    }
    stopifnot(is.numeric(dims), length(dims) == 2, 
        dims %in% seq_len(ncol(reducedDim(x, dr))))
        
    if (!all(color_by %in% rownames(x))) {
        stopifnot(length(color_by) == 1)
        if (!color_by %in% names(colData(x))) {
            .check_sce(x, TRUE)
            .check_pal(k_pal)
            .check_k(x, color_by)
            kids <- cluster_ids(x, color_by)
            nk <- nlevels(kids)
            if (length(k_pal) < nk)
                k_pal <- colorRampPalette(k_pal)(nk)
        } else kids <- NULL
    }
    
    # construct data.frame of reduced dimensions & relevant cell metadata
    xy <- reducedDim(x, dr)[, dims]
    colnames(xy) <- c("x", "y")
    df <- data.frame(colData(x), xy, check.names = FALSE)
    if (all(color_by %in% rownames(x))) {
        es <- as.matrix(assay(x, assay))
        es <- es[color_by, , drop = FALSE]
        if (scale) 
            es <- .scale_exprs(es, 1, q)
        df <- melt(
            cbind(df, t(es)), 
            id.vars = colnames(df))
        l <- switch(assay, exprs = "expression", assay)
        l <- paste0("scaled\n"[scale], l)
        scale <- scale_colour_gradientn(l, colors = a_pal)
        thm <- guide <- NULL
        color_by <- "value"
        facet <- facet_wrap("variable", ncol = ncol)
    } else if (is.numeric(df[[color_by]])) {
        if (scale) {
            vs <- as.matrix(df[[color_by]])
            df[[color_by]] <- .scale_exprs(vs, 2, q)
        }
        l <- paste0("scaled\n"[scale], color_by)
        scale <- scale_colour_gradientn(l, colors = a_pal)
        color_by <- sprintf("`%s`", color_by)
        facet <- thm <- guide <- NULL
    } else {
        facet <- NULL
        if (!is.null(kids)) {
            df[[color_by]] <- kids
            scale <- scale_color_manual(values = k_pal)
        } else scale <- NULL
        n <- nlevels(droplevels(factor(df[[color_by]])))
        guide <- guides(col = guide_legend(
            ncol = ifelse(n > 12, 2, 1),
            override.aes = list(alpha = 1, size = 3))) 
        thm <- theme(legend.key.height = unit(0.8, "lines"))
    }
    
    # set axes equal for linear dimension reductions
    if (dr %in% c("PCA", "MDS")) {
        asp <- coord_equal()
    } else asp <- NULL
    
    # get axes labels
    if (dr == "PCA") {
        labs <- paste0("PC", dims)
    } else labs <- paste(dr, "dim.", dims)
        
    # remove cells for which no reduced dimensions are available
    df <- df[!(is.na(df$x) | is.na(df$y)), ]
    
    p <- ggplot(df, aes_string("x", "y", col = color_by)) +
        geom_point(size = 0.4, alpha = 0.8) + 
        labs(x = labs[1], y = labs[2]) +
        facet + scale + guide + asp + 
        theme_minimal() + thm + theme(
            panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            axis.text = element_text(color = "black"),
            aspect.ratio = if (is.null(asp)) 1 else NULL)
    
    if (is.null(facet_by)) 
        return(p)
    
    if (is.null(facet)) {
        p + facet_wrap(facet_by, ncol = ncol)
    } else {
        if (nlevels(df$variable) == 1) {
            p + facet_wrap(facet_by, ncol = ncol) + 
                ggtitle(levels(df$variable))
        } else {
            fs <- c("variable", facet_by)
            ns <- vapply(df[fs], nlevels, numeric(1))
            if (ns[2] > ns[1]) fs <- rev(fs)
            p + facet_grid(reformulate(fs[1], fs[2]))
        }
    }
}
