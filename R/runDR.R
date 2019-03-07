#' @rdname runDR
#' @title Perform dim. reduction on a \code{daFrame}
#' 
#' @description Wrapper function to perform dimensionality reduction methods
#'   on a \code{daFrame} object using \code{scater}.
#' 
#' @param x a \code{\link{daFrame}}.
#' @param dr character string specifying the dimensionaly reduction method.
#' @param rows_to_use numeric vector of row indices (cells) to use.
#'   If NULL, all cells will be used. If a single integer value N,
#'   (default 1000) a subset of N cells will be drawn from each sample.
#' @param cols_to_use character vector in \code{colnames(x)} or numeric
#'   vector of column indices to use for computing reduced dimensions.
#' @param overwrite logical. Whether to force overwriting 
#'   any existing dimension reduction of type \code{dr}.
#' @param ... additional parameters passed to the dim. reduction method
#'   (see \code{\link[scater]{runTSNE}}, \code{\link[scater]{runPCA}}, 
#'   \code{\link[scater]{runMDS}}, \code{\link[scater]{runUMAP}}, 
#'   and \code{\link[scater]{runDiffusionMap}}).
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' daf <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' daf <- cluster(daf)
#' 
#' # PCA on all cells
#' daf <- runDR(daf, "PCA")
#' plotDR(daf, "PCA", color_by = "patient_id")
#' 
#' # UMAP on 1000 random cells
#' daf <- runDR(daf, "UMAP", rows_to_use = sample(nrow(daf), 1e3))
#' plotDR(daf, "UMAP", color_by = "condition")
#' 
#' # PCA on 200 cells per sample
#' set.seed(1)
#' daf <- runDR(daf, "PCA", rows_to_use = 200, overwrite = TRUE)
#' plotDR(daf, "PCA", color_by = "meta5")
#' 
#' # re-using PCA for t-SNE will fail when using different cells
#' \dontrun{
#' daf <- runDR(daf, "TSNE", rows_to_use = 1:500, use_dimred = "PCA")}
#' 
#' # use same seed to assure the same subset of cells is sampled
#' set.seed(1)
#' daf <- runDR(daf, "TSNE", rows_to_use = 200, use_dimred = "PCA")
#' plotDR(daf, "TSNE", color_by = "meta12")
#' 
#' # subsetting can be done 2-ways
#' identical(
#'   reducedDim(runDR(daf, "MDS", 1:100), "MDS"),
#'   reducedDim(runDR(daf[1:100, ], "MDS"), "MDS"))
#' 
#' @return a \code{daFrame} with an additional entry titled "\code{dr}"
#'   in the \code{reducedDims} slot of the input \code{daFrame}. 
#'   
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} 
#'   
#' @importFrom methods is
#' @importFrom scater runDiffusionMap runMDS runPCA runTSNE runUMAP
#' @importFrom SingleCellExperiment int_elementMetadata 
#'   reducedDimNames reducedDim reducedDims reducedDims<- 

setMethod("runDR", 
    signature = signature(x = "daFrame"),
    definition = function(x, 
        dr = c("TSNE", "PCA", "MDS", "UMAP", "DiffusionMap"), 
        rows_to_use = 1000, cols_to_use = NULL, overwrite = FALSE, ...) {
        
        dr <- match.arg(dr)
        stopifnot(is.logical(overwrite), length(overwrite) == 1)
        if (!overwrite && dr %in% reducedDimNames(x))
            stop(paste("A dimension reduction of type", dQuote(dr), 
                "is already present.\n  Run with 'overwrite = TRUE'",
                "to force replacement."))
            
        if (is.null(rows_to_use)) {
            rows_to_use <- seq_len(nrow(x))
        } else if (is.numeric(rows_to_use)) {
            stopifnot(all.equal(as.integer(rows_to_use), rows_to_use))
            if (length(rows_to_use) == 1) {
                idx <- split(seq_len(nrow(x)), sample_ids(x))
                rows_to_use <- lapply(idx, function(i)
                    sample(i, min(rows_to_use, length(i))))
                rows_to_use <- unname(unlist(rows_to_use))
            } else {
                stopifnot(all(rows_to_use %in% seq_len(nrow(x))))
            }
        } else {
            stop("Invalid argument 'rows_to_use'. Should be NULL,", 
                " a single integer value, or a vector of row indices.")
        }

        stopifnot(is.null(cols_to_use) || 
            is.character(cols_to_use) && all(cols_to_use %in% colnames(x)) ||
            is.numeric(cols_to_use) && all(cols_to_use %in% seq_len(ncol(x))))
        
        # default cols_to_use = type_markers(x)
        if (is.null(cols_to_use)) {
            if (length(type_markers(x)) < 3)
                stop("Please specify either which 'cols_to_use' or",
                    "\n  at least 3 'type' markers in 'colData(",
                    deparse(substitute(x)), ")$marker_class'")
            cols_to_use <- type_markers(x)
        }
        
        fun <- get(paste0("run", dr))
        rows_to_use <- sort(rows_to_use)
        y <- .rotate_daf(x[rows_to_use, ])
        if (dr == "TSNE" && !is.null(list(...)$use_dimred)) {
            use_dimred <- list(...)$use_dimred
            idx <- int_elementMetadata(x)
            k <- sprintf("idx.%s", use_dimred)
            k <- grep(k, names(idx))
            if (!identical(rows_to_use, which(idx[[k]])))
                stop(paste("The DR 'use_dimred =", dQuote(use_dimred), "'",
                    "seems to have been computed\n  for different cells",
                    "than specified with 'rows_to_use'.\n  Please run", 
                    sprintf("'runDR(.., %s)'", dQuote(use_dimred)), 
                    "using the same set of cells."))
            reducedDims(y) <- reducedDims(y)[use_dimred]
        } else {
            reducedDims(y) <- NULL
        }
        
        # compute DR
        y <- fun(y, exprs_values = "exprs", feature_set = cols_to_use, ...)
        # store coordinates in reducedDims
        res <- reducedDim(y, dr)
        colnames(res) <- paste0(dr, seq_len(ncol(res)))
        x@reducedDims[[dr]] <- res
        # store cells used in int_elementMetadata
        l <- seq_len(nrow(x)) %in% rows_to_use
        k <- sprintf("idx.%s", dr)
        x@int_elementMetadata[[k]] <- l
        return(x)
    }
)
