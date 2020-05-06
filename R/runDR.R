#' @rdname runDR
#' @title Dimension reduction
#' 
#' @description Wrapper around dimension reduction methods available 
#' through \code{scater}, with optional subsampling of cells per each sample.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param dr character string specifying which dimension reduction to use.
#' @param cells single numeric specifying the maximal number of cells
#'   per sample to use for dimension reduction; NULL for all cells.
#' @param features a character vector specifying which 
#'   antigens to use for dimension reduction; valid values are
#'   \code{"type"/"state"} for \code{type/state_markers(x)} 
#'   if \code{rowData(x)$marker_class} have been specified; 
#'   a subset of \code{rownames(x)}; NULL to use all features.
#' @param assay character string specifying which assay data to use
#'   for dimension reduction; valid values are \code{assayNames(x)}.
#' @param ... optional arguments for dimension reduction; passed to 
#'   \code{\link[scater]{runUMAP}}, \code{\link[scater]{runTSNE}}, 
#'   \code{\link[scater]{runPCA}}, \code{\link[scater]{runMDS}}
#'   and \code{\link[scater]{runDiffusionMap}}, respecttively.
#'   See \code{?"scater-red-dim-args"} for details.
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
#' # construct SCE
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run UMAP on <= 200 cells per sample
#' sce <- runDR(sce, features = type_markers(sce), cells = 100)
#' 
#' @importFrom scater runUMAP runTSNE runPCA runMDS runDiffusionMap
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
#' @importFrom SummarizedExperiment assayNames
#' @export

runDR <- function(x, 
    dr = c("UMAP", "TSNE", "PCA", "MDS", "DiffusionMap"), 
    cells = NULL, features = "type", assay = "exprs", ...) {
    
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    dr <- match.arg(dr)
    .check_assay(x, assay)
    fs <- .get_features(x, features)
    
    if (is.null(cells)) {
        # use all cells
        cs <- TRUE 
    } else {
        if (is.null(x$sample_id))
            stop("colData column sample_id not found,\n ", 
                " but is required to downsample cells.")
        stopifnot(
            is.numeric(cells), length(cells) == 1,
            as.integer(cells) == cells, cells > 0)
        # split cell indices by sample
        cs <- split(seq_len(ncol(x)), x$sample_id)
        # sample at most 'n' cells per sample
        cs <- unlist(lapply(cs, function(u)
            sample(u, min(cells, length(u)))))
    }
    
    # run dimension reduction
    fun <- get(paste0("run", dr))
    y <- fun(x[, cs], subset_row = fs, exprs_values = assay, ...)

    # return SCE when no cell subsetting has been done
    if (is.null(cells)) return(y)
    
    # else, write coordinates into original SCE
    xy <- reducedDim(y, dr)
    m <- matrix(NA, nrow = ncol(x), ncol = ncol(xy))
    m[cs, ] <- xy
    reducedDim(x, dr) <- m
    return(x)
}