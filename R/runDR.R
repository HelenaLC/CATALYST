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
#' 
#' # UMAP on 1000 random cells
#' daf <- runDR(daf, "UMAP", rows_to_use = sample(nrow(daf), 1e3))
#' 
#' reducedDims(daf)
#' head(reducedDim(daf, "UMAP"))
#' 
#' # PCA on 200 cells per sample
#' set.seed(1)
#' daf <- runDR(daf, "PCA", rows_to_use = 200, overwrite = TRUE)
#'
#' # re-using PCA for t-SNE will fail when using different cells
#' \dontrun{
#' daf <- runDR(daf, "TSNE", rows_to_use = 1:500, use_dimred = "PCA")}
#' 
#' # use same seed to assure the same subset of cells is sampled
#' set.seed(1)
#' daf <- runDR(daf, "TSNE", rows_to_use = 200, use_dimred = "PCA")
#' 
#' # number of rows used for each DR:
#' vapply(reducedDims(daf), nrow, numeric(1))
#' 
#' # running on subset can be done 2-ways
#' daf2 <- runDR(daf, "MDS", 1:100)
#' daf3 <- runDR(daf[1:100, ], "MDS")
#' 
#' # option 1 keeps object-dimensions
#' identical(dim(daf2), dim(daf))
#' 
#' # option 2 keeps only specified rows
#' all.equal(dim(daf3), c(100, ncol(daf)))
#' 
#' # reduced dimension are identical
#' identical(reducedDim(daf2), reducedDim(daf2))
#' 
#' # run t-SNE on B-cell clusters only
#' data(merging_table)
#' daf <- mergeClusters(daf, "meta20", merging_table, "merging")
#' cells_use <- grep("B-cells", cluster_ids(daf, "merging"))
#' daf <- runDR(daf, "TSNE", cells_use, overwrite = TRUE)
#' plotDR(daf, "TSNE", "merging")
#' 
#' @return a \code{daFrame} with an additional entry titled "\code{dr}"
#'   in the \code{reducedDims} slot of the input \code{daFrame}. 
#'   
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} 
#'   
#' @importFrom methods is
#' @importFrom scater runDiffusionMap runMDS runPCA runTSNE runUMAP
#' @importFrom SingleCellExperiment 
#'   int_elementMetadata int_elementMetadata<-
#'   reducedDim reducedDim<- reducedDims reducedDims<- 

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
        
        use_dimred <- list(...)$use_dimred
        stopifnot(use_dimred %in% reducedDimNames(x))
           
        if (is.null(rows_to_use)) {
            rows_to_use <- TRUE
        } else if (is.numeric(rows_to_use) && length(rows_to_use) == 1) {
            idx <- split(seq_len(nrow(x)), sample_ids(x))
            rows_to_use <- lapply(idx, function(i)
                sample(i, min(rows_to_use, length(i))))
            rows_to_use <- unname(unlist(rows_to_use))
        }
        is <- seq_len(nrow(x))[rows_to_use]
        rows_to_use <- sort(is)
        if (length(rows_to_use) == 0)
            stop("Argument 'rows_to_use' invalid; 0 cells specified.")

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
        y <- .rotate_daf(x[rows_to_use, ])
        if (is.null(use_dimred)) {
            reducedDims(y) <- NULL
        } else {
            if (dr != "PCA") {
                idx <- int_elementMetadata(x)
                k <- sprintf("idx.%s", use_dimred)
                k <- grep(k, names(idx))
                if (!identical(rows_to_use, which(idx[[k]])))
                    stop(paste("The DR 'use_dimred =", dQuote(use_dimred), "'",
                        "seems to have been computed\n  for different cells",
                        "than specified with 'rows_to_use'.\n  Please run",
                        sprintf("'runDR(.., %s)'", dQuote(use_dimred)),
                        "using the same set of cells."))
                i <- idx[[k]][rows_to_use]
                use_dr <- reducedDim(y, use_dimred)[i, ]
                use_dr <- as(list(use_dr), "SimpleList")
                names(use_dr) <- use_dimred
                reducedDims(y) <- use_dr
            } else if (dr == "PCA") {
                warning("'use_dimred' is not an argument to",
                    "  'scater::runPCA' and will be ignored.")
            }
        }
        
        # compute DR
        y <- fun(y, exprs_values = "exprs", feature_set = cols_to_use, ...)
        # store coordinates in reducedDims
        res <- reducedDim(y, dr)
        colnames(res) <- paste0(dr, seq_len(ncol(res)))
        x@reducedDims[[dr]] <- res
        # store cells used in int_elementMetadata
        l <- logical(nrow(x))
        l[rows_to_use] <- TRUE
        k <- sprintf("idx.%s", dr)
        int_elementMetadata(x)[[k]] <- l
        return(x)
    }
)
