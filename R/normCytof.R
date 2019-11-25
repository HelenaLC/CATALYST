#' @rdname normCytof
#' @title Bead-based normalization
#' 
#' @description 
#' an implementation of Finck et al.'s normalization of mass cytometry data 
#' using bead standards with automated bead gating.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param beads \code{"dvs"} (for bead masses 140, 151, 153 ,165, 175) or 
#'   \code{"beta"} (for bead masses 139, 141, 159, 169, 175) or 
#'   a numeric vector of masses.
#' @param out_path character string. If specified, outputs will be generated here. 
#'   If NULL (the default), \code{normCytof} will return a \code{SingleCellExperiment} 
#'   with an additional \code{assay} slot \code{"normed"} containing normalized data.
#'   If \code{remove_beads = TRUE}, bead events will be removed and returned separately.
#' @param fn
#'   a character string to use as the output file name. Defaults to 
#'   the file name of the input FCS file or \code{flowFrame}, respectively. 
#' @param fn_sep
#'   a character string to use to separate the output file name's prefix
#'   from the appendage.
#' @param remove_beads logical. If TRUE, bead events will be removed from
#'   the input \code{SingleCellExperiment} and returned as a separate object?
#' @param norm_to 
#'   a \code{\link{flowFrame}} or character of an FCS file from which baseline 
#'   values should be computed and to which the input data should be normalized.
#' @param k integer width of the median window used for bead smoothing.
#' @param trim a single non-negative numeric. 
#'   A \emph{median}+/-\code{trim}*\emph{mad} rule is applied to 
#'   preliminary bead populations to remove bead-bead doublets and 
#'   low signal beads prior to estimating normalization factors.
#' @param verbose logical. Should extra information on progress be reported?
#' @param plot logical. Should bead vs. DNA scatters and smoothed 
#'   bead intensities before vs. after normalization be plotted?
#' 
#' @return a list containing with the following elements: 
#' \itemize{
#' \item{a \code{SingleCellExperiment} with an additional assay slot \code{normed}}
#' \item{a \code{ggplot} of DNA vs. bead scatter plots}
#' \item{a \code{ggplot} of smoothed bead intensities before and after normalization}
#' }
#' 
#' @author Helena L. Crowell
#'
#' @references Finck, R. et al. (2013). 
#' Normalization of mass cytometry data with bead standards.
#' \emph{Cytometry A} \bold{83A}, 483-494.
#' 
#' @examples
#' data(raw_data)
#' sce <- fcs2sce(raw_data)
#' sce <- transform(sce)
#' res <- normCytof(sce, beads = "dvs", k = 80)
#' 
#' @importFrom flowCore colnames exprs isFCSfile read.FCS
#' @importFrom methods is
#' @importFrom Matrix colMeans
#' @importFrom matrixStats colAnys rowMins rowMedians rowMads
#' @importFrom stats approx runmed
#' @importFrom SingleCellExperiment int_metadata int_metadata<-
#' @importFrom SummarizedExperiment assayNames
#' @export

normCytof <- function(x, beads, assay = "exprs", k = 500, trim = 5, 
    remove_beads = TRUE, norm_to = NULL, out_path = NULL, verbose = TRUE) {
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        is.numeric(k), length(k) == 1, k > 1,
        is.null(out_path) || is.character(out_path) & 
            length(out_path) == 1 & dir.exists(out_path))
    
    # assure width of median window is odd
    if (k %% 2 == 0) k <- k + 1
    es <- assay(x, assay)

    # get time, length, DNA & bead channels
    ms <- .get_ms_from_chs(rownames(x))
    chs <- c(t = "time", l = "length", n = "filenum", dna = "Ir191|Ir193")
    chs <- sapply(chs, grep, rownames(x), ignore.case = TRUE, value = TRUE)
    n_beads <- length(chs$beads <- rownames(x)[
        .get_bead_cols(rownames(x), beads)])
    
    # identify bead singlets
    if (verbose) message("Identifying beads...")
    m <- match(unlist(chs[c("dna", "beads")]), rownames(x))
    key <- matrix(c(0, 0, rep(1, n_beads)), ncol = n_beads + 2, 
        dimnames = list("is_bead", ms[m]))
    key <- data.frame(key, check.names = FALSE)
    is_bead <- .get_bead_inds(x, key)
    
    # get all events that should be removed later
    # including bead-bead and cell-bead doublets
    ths <- rowMins(es[chs$beads, is_bead])
    rmv <- colSums(es[chs$beads, ] > ths) == n_beads
    x$remove <- rmv

    # trim tails
    z <- es[chs$beads, is_bead]
    meds <- rowMedians(z)
    mads <- rowMads(z) * trim
    diff <- abs(z - meds)
    ex <- colAnys(diff > mads)
    is_bead[which(is_bead)[ex]] <- FALSE
    x$is_bead <- is_bead

    if (verbose) message("Computing normalization factors...")

    # smooth bead intensitites by conversion to local medians
    smooth0 <- t(apply(es[chs$beads, is_bead], 1, runmed, k, "constant"))
    
    # get baseline (global mean)
    if (is.null(norm_to)) {
        bl <- rowMeans(es[chs$beads, is_bead])
    } else {
        if (is.character(norm_to)) {
            stopifnot(length(norm_to) == 1, isFCSfile(norm_to))
            ref <- read.FCS(norm_to, 
                transformation = FALSE, 
                truncate_max_range = FALSE)
        }
        bl <- colMeans(exprs(ref)[, .get_bead_cols(colnames(ref), y)])
    }
    
    # compute slopes (baseline versus smoothed bead intensitites)
    # & linearly interpolate slopes at non-bead events
    slopes <- colSums(smooth0 * bl) / colSums(smooth0 ^ 2)
    slopes <- approx(es[chs$t, is_bead], slopes, es[chs$t, ], rule = 2)$y

    # normalize raw bead intensities via multiplication with slopes
    i <- unlist(chs[c("t", "l", "n")])
    normed <-  es * slopes
    normed[i, ] <- es[i, ]
    assay(x, "normed") <- normed
   
    # smooth normalized beads
    smooth <- t(apply(normed[chs$beads, is_bead], 1, runmed, k, "constant"))

    ts <- es[chs$t, is_bead]
    ps <- list(
        scatter = .plot_bead_scatter(x, chs), 
        lines = .plot_smooth_beads(smooth0, smooth, ts))
    if (remove_beads) {
        z <- list(
            data = x[, !(is_bead | rmv)], 
            beads = x[, is_bead], 
            removed = x[, is_bead | rmv])
    } else {
        z <- list(x)
    }
    return(c(z, ps))
}