#' @rdname plotDR
#' @title Plot reduced dimensions
#' 
#' @description Dimension reduction plot colored 
#' by expression, cluster, sample or group ID.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param dr character string specifying which dimension reduction to use. 
#'   Should be one of \code{reducedDimNames(x)}; default to the 1st available.
#' @param color_by character string corresponding to a
#'   \code{colData(x)} column. Specifies the color coding.
#' @param facet_by character string specifying a
#'   non-numeric cell metadata column to facet by; 
#'   valid values are \code{names(colData(x))}.
#' @param k_pal character string specifying the cluster color palette; 
#'   ignored when \code{color_by} is not one of \code{names(cluster_codes(x))}. 
#'   If less than \code{nlevels(cluster_ids(x, k))} are supplied, colors will 
#'   be interpolated via \code{\link[grDevices:colorRamp]{colorRampPalette}}.
#' @param scale logical specifying whether expression should be scaled
#'   between 0 and 1 using lower (1\%) and upper (99\%) expression quantiles;
#'   ignored if \code{!all(color_by \%in\% rownames(x))}.
#' @param q single numeric in [0,0.5) determining the 
#'   quantiles to trim when \code{scale = TRUE}.
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
#' @importFrom grDevices colorRampPalette
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom SummarizedExperiment assay colData
#' @importFrom stats reformulate
#' @export

plotDR <- function(x, dr = NULL, 
    color_by = "condition", facet_by = NULL,
    k_pal = CATALYST:::.cluster_cols, scale = TRUE, q = 0.01) {
    
    # check validity of input arguments
    stopifnot(
        is(x, "SingleCellExperiment"),
        length(reducedDims(x)) != 0,
        is.logical(scale), length(scale) == 1,
        is.numeric(q), length(q) == 1, q >= 0, q < 0.5)
    .check_cd_factor(x, facet_by)
    
    if (!all(color_by %in% rownames(x))) {
        stopifnot(length(color_by) == 1)
        if (color_by %in% names(colData(x))) {
            .check_cd_factor(x, color_by)
        } else {
            .check_sce(x, TRUE)
            .check_pal(k_pal)
            .check_k(x, color_by)
            kids <- cluster_ids(x, color_by)
            nk <- nlevels(kids)
            if (length(k_pal) < nk)
                k_pal <- colorRampPalette(k_pal)(nk)
        }
    }

    if (is.null(dr)) {
        dr <- reducedDimNames(x)[1]
    } else {
        stopifnot(
            is.character(dr), length(dr) == 1, 
            dr %in% reducedDimNames(x))
    }
    
    # construct data.frame of reduced dimensions & relevant cell metadata
    df <- data.frame(colData(x), reducedDim(x, dr))
    if (all(color_by %in% rownames(x))) {
        es <- assay(x, "exprs")
        es <- es[color_by, , drop = FALSE]
        if (scale) 
            es <- .scale_exprs(es, 1, q)
        df <- melt(
            cbind(df, t(es)), 
            id.vars = colnames(df))
        scale <- scale_color_viridis_c(
            paste0("scaled\n"[scale], "expression"))
        thm <- guide <- NULL
        color_by <- "value"
        facet <- facet_wrap("variable")
        
    } else {
        facet <- NULL
        if (exists("kids")) {
            df[[color_by]] <- kids
            scale <- scale_color_manual(values = k_pal)
        } else scale <- NULL
        guide <- guides(col = guide_legend(
            override.aes = list(alpha = 1, size = 3))) 
        thm <- theme(legend.key.height  =  unit(0.8, "lines"))
    }
    
    # set axes equal for linear dimension reductions
    if (dr %in% c("PCA", "MDS")) {
        asp <- coord_equal()
    } else asp <- theme(aspect.ratio = 1)
    
    # get axes labels
    if (dr == "PCA") {
        labs <- paste(c("1st", "2nd"), "PC")
    } else labs <- paste(dr, "dim.", c(1, 2))
        
    # remove cells for which no reduced dimensions are available
    df <- df[!(is.na(df$X1) | is.na(df$X2)), ]
    
    p <- ggplot(df, aes_string("X1", "X2", col = color_by)) +
        geom_point(size = 0.4, alpha = 0.8) + 
        labs(x = labs[1], y = labs[2]) +
        facet + scale + guide + asp + 
        theme_minimal() + thm + theme(
            panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            axis.text = element_text(color = "black"))
    
    if (is.null(facet_by)) 
        return(p)
    
    if (is.null(facet)) {
        p + facet_wrap(facet_by)
    } else {
        p + facet_grid(reformulate("variable", facet_by))
    }
}
