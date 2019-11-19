#' @rdname assignPrelim
#' @title Single-cell debarcoding (1)
#' @description Assigns a preliminary barcode ID to each event.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param bc_key the debarcoding scheme. A binary matrix with sample names as row
#'   names and numeric masses as column names OR a vector of numeric masses 
#'   corresponding to barcode channels. When the latter is supplied, 
#'   `assignPrelim` will create a scheme of the appropriate format internally.
#' @param cf numeric. Cofactor used for asinh transformation.
#' @param verbose logical. Should extra information on progress be reported?
#' 
#' @return a \code{SingleCellExperiment} structured as follows: 
#' \describe{
#' \item{assays}{
#' \itemize{
#' \item{\code{counts} - raw counts}
#' \item{\code{exprs} - arcsinh-transformed counts}
#' \item{\code{scaled} - population-wise scaled 
#' expression using (95\%)-quantiles as boundaries}
#' }
#' }
#' \item{\code{colData}}{
#' \itemize{
#' \item{bc_id}{numeric verctor of barcode assignments}
#' \item{delta}{separation between positive and negative barcode populations}
#' }
#' }
#' \item{\code{metadata}}{
#' \itemize{
#' \item{bc_key}{the input debarcoding scheme}
#' }
#' }
#' }
#' 
#' @examples
#' data(sample_ff, sample_key)
#' es <- as.matrix(exprs(sample_ff))
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(
#'     assays = list(counts = t(es)),
#'     rowData = pData(parameters(sample_ff)))
#' sce <- assignPrelim(x = sce, y = sample_key)
#' 
#' @author Helena L. Crowell
#' 
#' @references 
#' Zunder, E.R. et al. (2015).
#' Palladium-based mass tag cell barcoding with a doublet-filtering scheme 
#' and single-cell deconvolution algorithm.
#' \emph{Nature Protocols} \bold{10}, 316-333. 
#' 
#' @importFrom matrixStats rowMaxs
#' @importFrom methods is
#' @importFrom stats quantile
#' @importFrom SingleCellExperiment altExp<-
#' @importFrom SummarizedExperiment assay assayNames
#' @export

assignPrelim <- function(x, bc_key, cf = 10, verbose = TRUE) {
    # check validity of input agruments
    stopifnot(
        is(x, "SingleCellExperiment"), 
        "counts" %in% assayNames(x),
        is.numeric(unlist(bc_key)), 
        is.vector(bc_key) | all(unlist(bc_key) %in% c(0, 1)),
        is.numeric(cf), length(cf) == 1, cf > 0,
        is.logical(verbose), length(verbose) == 1)
    
    if (is.vector(bc_key)) {
        n <- length(bc_key)
        bc_key <- matrix(diag(n), ncol = n, dimnames = list(bc_key, bc_key))
        bc_key <- data.frame(bc_key, check.names=FALSE)
    }
    
    # extract masses & check validity of debarcoding scheme
    n_bcs <- nrow(bc_key)
    ids <- rownames(bc_key)
    bc_ms <- as.numeric(colnames(bc_key))
    ms <- .get_ms_from_chs(rownames(x))
    if (any(!bc_ms %in% ms))
        stop("Couldn't match masses extracted from", 
            " channel names and debarcoding scheme.")
    
    # get columns corresponding to barcode channels
    bc_cols <- vapply(bc_ms, function(u) which(ms == u), numeric(1))
    if (length(bc_cols) != ncol(bc_key))
        stop("Not all barcode channels found.")
    
    # subset & transform barcode channels
    assay(x, "exprs") <- asinh(assay(x, "counts") / cf)
    bc_es <- assay((y <-  x[bc_cols, ]), "exprs")
     
    # assign barcode ID to each cell
    if (verbose) message("Debarcoding data...")
    bc_ids <- .get_ids(bc_es, bc_key, ids, verbose)
    
    # rescale transformed barcodes for each population 
    # using preliminary assignments
    if (verbose) message("Normalizing...")
    normed_bcs <- matrix(0, 
        nrow=nrow(bc_es), ncol=ncol(x), 
        dimnames=list(rownames(bc_es), NULL))
    pos <- lapply(ids, `==`, bc_ids)
    for (i in seq_along(ids))
        if (any(pos[[i]])) {
            pos_bcs <- bc_es[bc_key[i, ] == 1, pos[[i]]]
            q95 <- quantile(pos_bcs, 0.95)
            normed_bcs[, pos[[i]]] <- bc_es[, pos[[i]]]/q95
        }
    
    # get deltas from normalized intensities 
    if (verbose) message("Computing deltas...")
    deltas <- .get_deltas(normed_bcs, bc_key, verbose)
    
    assay(y, "scaled") <- normed_bcs
    y$bc_id <- bc_ids
    y$delta <- deltas
    metadata(y)$bc_key <- bc_key
    altExp(x, "barcodes") <- y
    return(x)
}