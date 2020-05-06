#' @rdname assignPrelim
#' @title Single-cell debarcoding (1)
#' @description Assigns a preliminary barcode ID to each event.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param bc_key the debarcoding scheme. A binary matrix with sample names as
#'   row names and numeric masses as column names OR a vector of numeric masses 
#'   corresponding to barcode channels. When the latter is supplied, 
#'   `assignPrelim` will create a scheme of the appropriate format internally.
#' @param assay character string specifying which assay to use.
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
#' sce <- prepData(sample_ff)
#' sce <- assignPrelim(sce, sample_key)
#' table(sce$bc_id)
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
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

assignPrelim <- function(x, bc_key, assay = "exprs", verbose = TRUE) {
    # check validity of input arguments
    args <- as.list(environment())
    .check_args_assignPrelim(args)
    
    if (is.vector(bc_key)) {
        n <- length(bc_key)
        bc_key <- matrix(diag(n), ncol = n, dimnames = list(bc_key, bc_key))
        bc_key <- data.frame(bc_key, check.names = FALSE)
    }
    
    # extract masses & check validity of debarcoding scheme
    n_bcs <- length(ids <- rownames(bc_key))
    bc_ms <- as.numeric(colnames(bc_key))
    if (any(is.na(bc_ms)))
        stop("Column names of 'bc_key' should be numeric.")
    chs <- channels(x)
    ms <- .get_ms_from_chs(chs)
    if (any(!bc_ms %in% ms))
        stop("Couldn't match masses extracted from", 
            " channel names and debarcoding scheme.")
    
    # get columns corresponding to barcode channels
    bc_chs <- vapply(bc_ms, function(u) match(u, ms), numeric(1))
    if (length(bc_chs) != ncol(bc_key))
        stop("Not all barcode channels found.")
    
    # specify barcode channels
    rowData(x)$is_bc <- is_bc <- seq_len(nrow(x)) %in% bc_chs
    
    # assign barcode ID to each cell
    if (verbose) message("Debarcoding data...")
    x$bc_id <- .get_ids(assay(x, assay)[bc_chs, ], bc_key, ids, verbose)
    
    # rescale transformed barcodes for each population 
    # using preliminary assignments
    if (verbose) message("Normalizing...")
    # split cells by populations
    cs <- split(seq_len(ncol(x)), x$bc_id)[ids]
    scaled <- matrix(NA, nrow(x), ncol(x), dimnames = dimnames(x))
    scaled[is_bc, unlist(cs[ids])] <- do.call("cbind", 
        lapply(ids, function(id) {
            y <- assay(x, assay)[, cs[[id]], drop = FALSE]
            pos <- y[bc_chs[bc_key[id, ] == 1], ]
            y[is_bc, ] / quantile(pos, 0.95)
        })
    )
    assay(x, "scaled", withDimnames = FALSE) <- scaled
    
    # get deltas from normalized intensities 
    if (verbose) message("Computing deltas...")
    y <- assay(x, "scaled", withDimnames = FALSE)[is_bc, ]
    x$delta <- .get_deltas(y, bc_key, verbose)
    
    # store debarcoding scheme in metadata & return SCE
    metadata(x)$bc_key <- bc_key
    return(x)
}
