#' @rdname tSNE
#' @title Run t-SNE
#' 
#' @description Runs t-SNE dimensionality reduction on a \code{\link{daFrame}}.
#' 
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param n 
#'   numeric. Specifies the number of cells to downsample to per sample.
#' @param seed 
#'   numeric. Specifies the seed to be set before sampling
#' 
#' @return Writes the tSNE coordinates and the indicies of the events used 
#' for their computation into the metadata slot of the input \code{daFrame}.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
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
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- tSNE(re, n=250)
#'
#' @import Rtsne SummarizedExperiment
# ------------------------------------------------------------------------------

setMethod(f="tSNE",
    signature=signature(x="daFrame"),
    definition=function(x, n=1000, seed=42) {
        message("o downsampling to ", n, " events per sample...")
        dups <- which(!duplicated(exprs(x)[, type_markers(x)]))
        n_cells <- pmin(metadata(x)$n_cells, n)
        inds <- split(seq_len(nrow(exprs(x))), sample_ids(x))
        set.seed(seed)
        tsne_inds <- lapply(names(inds), function(i) {
            s <- sample(inds[[i]], n_cells[i], replace=FALSE)
            intersect(s, dups)
        })
        message("o running tSNE...")
        tsne_inds <- unlist(tsne_inds)
        tsne_es <- exprs(x)[tsne_inds, cols]
        tsne <- Rtsne(tsne_es, check_duplicates=FALSE, pca=FALSE)
        metadata(x)$tsne <- tsne
        metadata(x)$tsne_inds <- tsne_inds
        return(x)
    }
)
