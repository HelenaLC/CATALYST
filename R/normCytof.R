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
#' @param assay character string specifying which assay data to use.
#'   Should be one of \code{assayNames(x)}.
#' @param remove_beads logical. If TRUE, bead events will be removed from
#'   the input \code{SingleCellExperiment} and returned as a separate object?
#' @param norm_to a \code{\link[flowCore]{flowFrame}} or path to
#'   an FCS file from which to compute baseline bead intensities,
#'   and to which the input data should be normalized to.
#' @param k integer width of the median window used 
#'   for bead smoothing (affects visualizations only!).
#' @param trim a single non-negative numeric. 
#'   A \emph{median}+/-\code{trim}*\emph{mad} rule is applied to 
#'   preliminary bead populations to remove bead-bead doublets and 
#'   low signal beads prior to estimating normalization factors.
#' @param verbose logical. Should extra information on progress be reported?
#' @param plot logical. Should bead vs. DNA scatters and smoothed bead 
#'   intensities before vs. after normalization be included in the output?
#' 
#' @return a list of the following \code{SingleCellExperiment}... 
#' \itemize{
#' \item{\code{data}: SCE with additional assay slot \code{normed} containing 
#'   normalized (and optionally filtered, if \code{remove_beads = TRUE}) data}
#' \item{\code{beads}, \code{removed}: SCEs containing subsets of events 
#'   identified as beads and that were removed, respectively. 
#'   The latter includes bead-cell and cell-cell doublets)}
#' } ...and \code{ggplot} objects:
#' \itemize{
#' \item{\code{scatter}: scatter plot of DNA vs. bead 
#'   intensities with indication of the applied gates}
#' \item{\code{lines}: running-median smoothed bead 
#'   intensities before and after normalization}
#' }
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#'
#' @references Finck, R. et al. (2013). 
#' Normalization of mass cytometry data with bead standards.
#' \emph{Cytometry A} \bold{83A}, 483-494.
#' 
#' @examples
#' data(raw_data)
#' sce <- fcs2sce(raw_data)
#' res <- normCytof(sce, beads = "dvs", k = 80)
#' 
#' ncol(res$beads)   # no. of bead events
#' ncol(res$removed) # no. of events removed
#' 
#' res$scatter # plot DNA vs. bead intensities including applied gates
#' res$lines   # plot smoothed bead intensities before vs. after normalization
#' 
#' @importFrom flowCore colnames exprs isFCSfile read.FCS
#' @importFrom methods is
#' @importFrom Matrix colMeans
#' @importFrom matrixStats colAnys rowMins rowMedians rowMads
#' @importFrom stats approx runmed
#' @importFrom SingleCellExperiment int_colData int_metadata int_metadata<-
#' @importFrom SummarizedExperiment assayNames assay assay<-
#' @export

normCytof <- function(x, beads, assay = "exprs", k = 500, trim = 5, 
    remove_beads = TRUE, norm_to = NULL, plot = TRUE, verbose = TRUE) {
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        is.numeric(k), length(k) == 1, k > 1,
        is.numeric(trim), length(trim) == 1, trim >= 0,
        is.logical(remove_beads), length(remove_beads) == 1,
        is.null(norm_to) || is(norm_to, "flowFrame") 
        || is.character(norm_to) & length(norm_to) == 1 & isFCSfile(norm_to),
        is.logical(plot), length(plot) == 1,
        is.logical(verbose), length(verbose) == 1,
        any(grepl("Ir191|Ir193", (chs <- rownames(x)), ignore.case = TRUE)),
        any(grepl("time", names(icd <- int_colData(x)), ignore.case = TRUE)))

    # get times, DNA & bead channels
    ts <- icd[[grep("time", names(icd), ignore.case = TRUE)]]
    dna_chs <- grep("Ir191|Ir193", chs, ignore.case = TRUE, value = TRUE)
    bead_chs <- chs[.get_bead_cols(chs, beads)]
    rowData(x)$bead_ch <- chs %in% bead_chs
    n_beads <- length(bead_chs)
    
    # identify bead singlets
    if (verbose) message("Identifying beads...")
    ms <- .get_ms_from_chs(chs)
    m <- match(c(dna_chs, bead_chs), chs)
    key <- matrix(c(0, 0, rep(1, n_beads)), ncol = n_beads + 2, 
        dimnames = list("is_bead", ms[m]))
    key <- data.frame(key, check.names = FALSE)
    is_bead <- .get_bead_inds(x, key, assay = assay)
    
    # subset specified assay
    es <- assay(x, assay)
    
    # get all events that should be removed later
    # including bead-bead and cell-bead doublets
    ths <- rowMins(es[bead_chs, is_bead])
    rmv <- colSums(es[bead_chs, ] > ths) == n_beads
    x$remove <- rmv

    # trim tails
    z <- es[bead_chs, is_bead]
    meds <- rowMedians(z)
    mads <- rowMads(z) * trim
    diff <- abs(z - meds)
    ex <- colAnys(diff > mads)
    is_bead[which(is_bead)[ex]] <- FALSE
    x$is_bead <- is_bead
    
    # get baselines (global mean)
    if (is.null(norm_to)) {
        bl <- rowMeans(es[bead_chs, is_bead])
    } else {
        if (is.character(norm_to)) {
            ref <- read.FCS(norm_to, 
                transformation = FALSE, 
                truncate_max_range = FALSE)
        } else {
            ref <- norm_to
        }
        bl <- colMeans(exprs(ref)[, .get_bead_cols(colnames(ref), beads)])
    }
    
    # smooth bead intensitites by conversion to local medians
    if (verbose) message("Computing normalization factors...")
    
    # assure width of median window is odd 
    if (k %% 2 == 0) k <- k + 1
    smooth0 <- t(apply(es[bead_chs, is_bead], 1, runmed, k, "constant"))

    # compute slopes (baseline versus smoothed bead intensitites)
    # & linearly interpolate slopes at non-bead events
    slopes <- colSums(smooth0 * bl) / colSums(smooth0 ^ 2)
    slopes <- approx(ts[is_bead], slopes, ts, rule = 2)$y
    
    # normalize raw bead intensities via multiplication with slopes
    assay(x, "normed", withDimnames = FALSE) <- sweep(es, 2, slopes, "*")
   
    # smooth normalized beads
    y <- assay(x, "normed", withDimnames = FALSE)[bead_chs, is_bead]
    smooth <- t(apply(y, 1, runmed, k, "constant"))

    ps <- NULL
    if (plot) ps <- list(
        scatter = .plot_bead_scatter(x, dna_chs, bead_chs, assay), 
        lines = .plot_smooth_beads(smooth0, smooth, ts[is_bead]))
    if (remove_beads) {
        z <- list(
            data = x[, !(is_bead | rmv)], 
            beads = x[, is_bead], 
            removed = x[, is_bead | rmv])
    } else {
        z <- list(data = x)
    }
    return(c(z, ps))
}
