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
#' @param assay character string specifying which assay to use. Should 
#'   correspond to count-like data, as linearity-assumptions underlying 
#'   spillover estimation won't hold for non-linearly transformed data.
#' @param transform logical specifying whether compensated data 
#'   should be arcsinh-transformed in which case an additional
#'   assay \code{"exprs_comped"} will be added to the output SCE.
#' @param cofactor cofactor to use for optional arcsinh-transformation.
#'   If NULL, \code{compCytof} will try and access \code{metadata(x)$cofactor}.
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
#' sce <- fcs2sce(ss_exp, by_time = FALSE)
#' bc_ms <- c(139, 141:156, 158:176)
#' sce <- assignPrelim(x = sce, bc_key = bc_ms)
#' sce <- estCutoffs(x = sce)
#' sce <- applyCutoffs(x = sce)
#' sce <- computeSpillmat(x = sce)
#' (sce <- compCytof(x = sce))
#' 
#' library(SingleCellExperiment)
#' i <- "Dy162Di"; j <- "Dy163Di"
#' par(mfrow = c(1, 2))
#' for (a in c("exprs", "exprs_comped"))
#'   plot(main = a, xlab = i, ylab = j,
#'     assay(sce[i, ], a), assay(sce[j, ], a))
#'
#' @importFrom methods is
#' @importFrom nnls nnls
#' @importFrom S4Vectors metadata
#' @importFrom flowCore compensate exprs
#' @importFrom SummarizedExperiment assay assay<- assayNames
#' @export

compCytof <- function(x, sm = NULL, method = c("nnls", "flow"),
    assay = "counts", transform = TRUE, cofactor = NULL, 
    isotope_list = CATALYST::isotope_list) {
    
    # check validity of input arguments
    method <- match.arg(method)
    stopifnot(is(x, "SingleCellExperiment"),
        !is.null(sm) || !is.null(metadata(x)$spillover_matrix),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        !is.null(cofactor) | !is.null(cofactor <- metadata(x)$cofactor),
        is.numeric(cofactor), length(cofactor) == 1, cofactor > 0)
    
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
    a <- sprintf("%s.%scomped", assay, method)
    assay(x, a) <- y
    
    # (optionally) apply arcsinh-transformation to compensated data
    if (transform) {
        a <- sprintf("exprs.%scomped", method)
        assay(x, a) <- asinh(y/cofactor)
    }
    return(x)
}
