#' @rdname compCytof
#' @title Compensate CyTOF data
#' 
#' @description Compensates a mass spectrometry based experiment using a
#' provided spillover matrix & assuming a linear spillover in the experiment.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}} 
#'   OR a character string specifying the location of 
#'   FCS files that should be compensates.
#' @param sm a spillover matrix.
#' @param method \code{"flow"} or \code{"nnls"}.
#' @param assay character string specifying which assay data to use; 
#'   should be one of \code{assayNames(x)} and correspond to count-like data, 
#'   as linearity assumptions underlying compensation won't hold otherwise.
#' @param overwrite logical; should the specified \code{assay} slot 
#'   (and \code{exprs}, when \code{transform = TRUE}) be overwritten 
#'   with the compensated data? If \code{FALSE}, compensated counts 
#'   (and expressions, if \code{transform = TRUE}) will be stored in 
#'   assay(s) \code{compcounts/exprs}, respectively.
#' @param transform logical; should normalized counts be 
#'   arcsinh-transformed with the specified \code{cofactor}(s)?
#' @param cofactor numeric cofactor(s) to use for optional 
#'   arcsinh-transformation when \code{transform = TRUE};
#'   single value or a vector with channels as names.
#'   If NULL, \code{compCytof} will try and access the cofactor(s)
#'   stored in \code{int_metadata(x)}, thus re-using the same 
#'   transformation applied previously.
#' @param isotope_list named list. Used to validate the input spillover matrix.
#'   Names should be metals; list elements numeric vectors of their isotopes.
#'   See \code{\link{isotope_list}} for the list of isotopes used by default.
#' 
#' @details
#' If the spillover matrix (SM) does not contain the same set of columns as 
#' the input experiment, it will be adapted according to the following rules:
#' \enumerate{
#' \item{columns present in the SM but not in the input data 
#' will be removed from it}
#' \item{non-metal columns present in the input but not in the SM 
#' will be added such that they do neither receive nor cause spill}
#' \item{metal columns that have the same mass as a channel present in the SM 
#' will receive (but not emit) spillover according to that channel}
#' \item{if an added channel could potentially receive spillover (as it has 
#' +/-1M or +16M of, or is of the same metal type as another channel measured), 
#' a warning will be issued as there could be spillover interactions that have
#' been missed and may lead to faulty compensation}}
#' 
#' @return 
#' Compensates the input \code{\link{flowFrame}} or, 
#' if \code{x} is a character string, all FCS files in the specified location. 
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch} 
#' & Vito RT Zanotelli
#' 
#' @examples
#' # deconvolute single-stained control samples
#' data(ss_exp)
#' sce <- prepData(ss_exp)
#' bc_ms <- c(139, 141:156, 158:176)
#' sce <- assignPrelim(sce, bc_ms)
#' sce <- applyCutoffs(estCutoffs(sce))
#' 
#' # estimate spillover matrix 
#' sce <- computeSpillmat(sce)
#' 
#' # compensate & store compensated data in separate assays
#' sce <- compCytof(sce, overwrite = FALSE)
#' assayNames(sce)
#' 
#' # biscatter before vs. after compensation
#' chs <- c("Dy162Di", "Dy163Di")
#' m <- match(chs, channels(sce))
#' i <- rownames(sce)[m][1]
#' j <- rownames(sce)[m][2]
#' 
#' par(mfrow = c(1, 2))
#' for (a in c("exprs", "compexprs")) {
#'   es <- assay(sce, a)
#'   plot(es[i, ], es[j, ], cex = 0.2, pch = 19,
#'        main = a, xlab = i, ylab = j)
#' }
#'
#' @importFrom nnls nnls
#' @importFrom S4Vectors metadata
#' @importFrom flowCore compensate exprs
#' @importFrom SingleCellExperiment int_metadata
#' @importFrom SummarizedExperiment assay assay<- assayNames rowData
#' @export

compCytof <- function(x, sm = NULL, method = c("nnls", "flow"),
    assay = "counts", overwrite = TRUE, transform = TRUE, cofactor = NULL, 
    isotope_list = CATALYST::isotope_list) {
    # check validity of input arguments
    method <- match.arg(method)
    args <- as.list(environment())
    .check_args_compCytof(args)
    if (is.null(cofactor))
        cofactor <- int_metadata(x)$cofactor
    
    chs0 <- rownames(x)
    rownames(x) <- channels(x)
    if (is.null(sm)) sm <- metadata(x)$spillover_matrix
    if (!is.matrix(sm)) sm <- as.matrix(sm)
    suppressMessages(sm <- adaptSpillmat(sm, rownames(x), isotope_list))
    
    # apply compensation & store compensated data in assays
    y <- switch(method, 
        flow = {
            a <- as.matrix(assay(x, assay))
            ff <- flowFrame(t(a))
            ff <- compensate(ff, sm)
            t(exprs(ff))
        },
        nnls = apply(assay(x, assay), 2, 
            function(u) nnls(t(sm), u)$x))
    c <- ifelse(overwrite, assay, "compcounts")
    assay(x, c, FALSE) <- y
    
    # do arcsinh-transformation on compensated counts
    if (transform) {
        e <- ifelse(overwrite, "exprs", "compexprs")
        x <- .transform(x, cofactor, ain = c, aout = e)
    }
    rownames(x) <- chs0
    return(x)
}
