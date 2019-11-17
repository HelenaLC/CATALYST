#' @rdname assignPrelim
#' @title Single-cell debarcoding (1)
#' @description Assigns a preliminary barcode ID to each event.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param y the debarcoding scheme. A binary matrix with sample names as row
#'   names and numeric masses as column names OR a vector of numeric masses 
#'   corresponding to barcode channels. When the latter is supplied, 
#'   `assignPrelim` will create a scheme of the appropriate format internally.
#' @param cofactor numeric. Cofactor used for asinh transformation.
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
#' @importFrom SummarizedExperiment assay assayNames
#' @export

assignPrelim <- function(x, y, cofactor=10, verbose=TRUE) {
    # check validity of input agruments
    stopifnot(is(x, "SingleCellExperiment"), "counts" %in% assayNames(x),
        is.numeric(unlist(y)), all(unlist(y) %in% c(0, 1)),
        is.numeric(cofactor), length(cofactor) == 1, cofactor > 0,
        is.logical(verbose), length(verbose) == 1)
    
    if (is.vector(y)) {
        n <- length(y)
        y <- matrix(diag(n), ncol = n, dimnames = list(y, y))
        y <- data.frame(y, check.names=FALSE)
    }
    
    # extract masses & check validity of debarcoding scheme
    n_bcs <- nrow(y)
    ids <- rownames(y)
    bc_ms <- as.numeric(colnames(y))
    ms <- .get_ms_from_chs(rownames(x))
    if (any(!bc_ms %in% ms))
        stop("Couldn't match masses extracted from", 
            " channel names and debarcoding scheme.")
    
    # get columns corresponding to barcode channels
    bc_cols <- vapply(bc_ms, function(u) which(ms == u), numeric(1))
    if (length(bc_cols) != ncol(y))
        stop("Not all barcode channels found.")
    
    # subset & transform barcode channels
    x <- x[bc_cols, ]
    es <- assay(x, "exprs") <- asinh(assay(x, "counts") / cofactor)
     
    # assign barcode ID to each cell
    if (verbose) message("Debarcoding data...")
    bc_ids <- .get_ids(es, y, ids, verbose)
    
    # rescale transformed barcodes for each population 
    # using preliminary assignments
    if (verbose) message("Normalizing...")
    normed_bcs <- matrix(0, 
        nrow=nrow(es), ncol=ncol(x), 
        dimnames=list(rownames(es), NULL))
    pos <- lapply(ids, `==`, bc_ids)
    for (i in seq_along(ids))
        if (any(pos[[i]])) {
            pos_bcs <- es[y[i, ] == 1, pos[[i]]]
            q95 <- quantile(pos_bcs, 0.95)
            normed_bcs[, pos[[i]]] <- es[, pos[[i]]]/q95
        }
    assay(x, "scaled") <- normed_bcs
    
    # get deltas from normalized intensities 
    if (verbose) message("Computing deltas...")
    deltas <- .get_deltas(normed_bcs, y, verbose)
    
    # compute well-wise yield for each separation cutoff
    if (verbose) message("Computing counts and yields...")
    n_seps <- length(seps <- seq(0, 1, 0.01))
    yields <- counts <- matrix(0, nrow = n_bcs, ncol = n_seps)
    for (i in seq_along(ids)) {
        #pos <- which(inds == i)
        for (j in seq_along(seps)) {
            k <- deltas[pos[[i]]] >= seps[j]
            yields[i, j] <- sum(k)
            counts[i, j] <- sum(k & deltas[pos[[i]]] < seps[j + 1])
            if (j == n_seps)
                counts[i, j] <- sum(k)
        }
    }
    
    # normalize yields
    norm_val <- rowMaxs(yields)
    norm_val[norm_val == 0] <- 1
    yields <- yields / norm_val
    
    rownames(counts) <- rownames(yields) <- ids
    colnames(counts) <- colnames(yields) <- seps
    
    x$bc_id <- bc_ids
    x$delta <- deltas
    metadata(x)$bc_key <- y
    return(x)
}