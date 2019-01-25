#' @rdname tSNE
#' @title Run t-SNE
#' 
#' @description Runs t-SNE dimensionality reduction on a \code{\link{daFrame}}.
#' 
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param cols_to_use
#'   a character vector. Specifies which antigens to use for clustering.
#'   If NULL, the function will attempt to use the \code{type_markers(x)}.
#' @param n 
#'   numeric. Specifies the number of cells to downsample to per sample.
#' @param verbose
#'   logical. Should information on progress be reported?
#' @param seed 
#'   numeric. Specifies the seed to be set before sampling
#' 
#' @return Writes the tSNE coordinates and the indicies of the events used 
#' for their computation into the metadata slot of the input \code{daFrame}.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @seealso \code{\link{plotSNE}}
#' 
#' @examples
#' # construct daFrame
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run t-SNE
#' re <- tSNE(re, n=50)
#' 
#' par(pty="s")
#' tsne <- S4Vectors::metadata(re)$tsne$Y
#' plot(tsne, pch=20)
#'
#' @import Rtsne  
#' @importFrom S4Vectors metadata
# ------------------------------------------------------------------------------

setMethod(f="tSNE",
    signature=signature(x="daFrame"),
    definition=function(x, cols_to_use=NULL, n=1000, verbose=TRUE, seed=42) {
        if (is.null(cols_to_use)) {
            type_markers <- colData(x)$marker_class == "type"
            if (sum(colData(x)$marker_class == "type") < 3)
                stop("Please specify either which 'cols_to_use' or", 
                    "\nat least 3 'type' markers in 'colData(", 
                    deparse(substitute(x)), ")$marker_class'")
            cols_to_use <- type_markers
        }
        
        if (verbose)
            message("o downsampling to ", n, " events per sample...")
        dups <- which(!duplicated(exprs(x)[, cols_to_use]))
        n_cells <- pmin(metadata(x)$n_cells, n)
        inds <- split(seq_len(nrow(exprs(x))), sample_ids(x))
        set.seed(seed)
        tsne_inds <- lapply(names(inds), function(i) {
            s <- sample(inds[[i]], n_cells[i], replace=FALSE)
            intersect(s, dups)
        })
        if (verbose)
            message("o running tSNE...")
        tsne_inds <- unlist(tsne_inds)
        tsne_es <- exprs(x)[tsne_inds, cols_to_use]
        set.seed(seed)
        tsne <- Rtsne(tsne_es, check_duplicates=FALSE, pca=FALSE)
        metadata(x)$tsne <- tsne
        metadata(x)$tsne_inds <- tsne_inds
        return(x)
    }
)
