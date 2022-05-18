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
#' @param dna numeric vector of masses corresponding to DNA channels  
#'   (only one is required; output scatter plot (see Value section)  
#'   will be generated using the first matching channel).
#' @param assays lnegth 2 character string specifying 
#'   which assay data to use; both should be in \code{assayNames(x)} 
#'   and correspond to count- and expression-like data, respectively.
#' @param remove_beads logical. If TRUE, bead events will be removed from
#'   the input \code{SingleCellExperiment} and returned as a separate object?
#' @param norm_to 
#'   a \code{\link[flowCore:flowFrame-class]{flowFrame}} or character string
#'   specifying an FCS file from which to compute baseline bead intensities,
#'   and to which the input data should be normalized to.
#' @param k integer width of the median window used 
#'   for bead smoothing (affects visualizations only!).
#' @param trim a single non-negative numeric. 
#'   A \emph{median}+/-\code{trim}*\emph{mad} rule is applied to 
#'   preliminary bead populations to remove bead-bead doublets and 
#'   low signal beads prior to estimating normalization factors.
#' @param overwrite logical; should the specified \code{assays}
#'   (both, when \code{transform = TRUE}) be overwritten 
#'   with the normalized data? If \code{FALSE}, normalized counts 
#'   (and expressions, if \code{transform = TRUE}) will be stored in 
#'   assay(s) \code{normcounts/exprs}, respectively.
#' @param transform logical; should normalized counts be 
#'   arcsinh-transformed with the specified \code{cofactor}(s)?
#' @param cofactor numeric cofactor(s) to use for optional 
#'   arcsinh-transformation when \code{transform = TRUE};
#'   single value or a vector with channels as names.
#'   If NULL, \code{normCytof} will try and access the cofactor(s)
#'   stored in \code{int_metadata(x)}, thus re-using the same 
#'   transformation applied previsouly.
#' @param plot logical; should bead vs. DNA scatters and smoothed bead 
#'   intensities before vs. after normalization be included in the output?
#' @param verbose logical; should extra information on progress be reported?
#' 
#' @return a list of the following \code{SingleCellExperiment}... 
#' \itemize{
#' \item{\code{data}: 
#'   The filtered input SCE (when \code{remove_beads = TRUE}); 
#'   otherwise, \code{colData} columns \code{is_bead} and \code{remove} 
#'   indicate whether an event as been identified as a bead or doublet.
#'   If \code{overwrite = FALSE}, assays \code{normcounts/exprs} are added;
#'   otherwise, the specified \code{counts/exprs} assays are overwritten.}
#' \item{\code{beads}, \code{removed}: 
#'   SCEs containing subsets of events identified as beads 
#'   and that were removed, respectively. The latter includes 
#'   bead-cell and cell-cell doublets)}
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
#' sce <- prepData(raw_data)
#' 
#' # apply normalization & write normalized data to separate assays
#' res <- normCytof(sce, beads = "dvs", k = 80, overwrite = FALSE) 
#' 
#' ncol(res$beads)   # no. of bead events
#' ncol(res$removed) # no. of events removed
#' 
#' res$scatter # plot DNA vs. bead intensities including applied gates
#' res$lines   # plot smoothed bead intensities before vs. after normalization
#' 
#' # filtered SCE now additionally includes 
#' # normalized count & expression data
#' assayNames(res$data)
#' 
#' @importFrom flowCore colnames exprs read.FCS
#' @importFrom Matrix colMeans
#' @importFrom matrixStats colAnys rowMins rowMedians rowMads
#' @importFrom stats approx runmed
#' @importFrom SingleCellExperiment int_colData int_metadata
#' @importFrom SummarizedExperiment assayNames assay assay<- rowData
#' @export

normCytof <- function(x, 
    beads = c("dvs", "beta"), 
    dna = c(191, 193), k = 500, trim = 5, 
    remove_beads = TRUE, norm_to = NULL, 
    assays = c("counts", "exprs"),
    overwrite = TRUE, transform = TRUE, cofactor = NULL, 
    plot = TRUE, verbose = TRUE) {
    # check validity of input arguments
    args <- as.list(environment())
    .check_args_normCytof(args)
    if (is.null(cofactor))
        cofactor <- int_metadata(x)$cofactor
    if (is.character(beads))
        beads <- match.arg(beads)

    # get times, DNA & bead channels
    icd <- int_colData(x)
    chs <- channels(x)
    chs0 <- rownames(x)
    rownames(x) <- chs
    stopifnot(
        any(grepl("Ir191|Ir193", chs, ignore.case = TRUE)),
        any(grepl("time", names(icd), ignore.case = TRUE)))
    
    ts <- icd[[grep("time", names(icd), ignore.case = TRUE)]]
    n_dna <- length(dna_chs <- chs[.get_dna_cols(chs, dna)])
    n_beads <- length(bead_chs <- chs[.get_bead_cols(chs, beads)])
    rowData(x)$bead_ch <- chs %in% bead_chs
    
    # identify bead singlets
    if (verbose) message("Identifying beads...")
    ms <- .get_ms_from_chs(chs)
    m <- match(c(dna_chs, bead_chs), chs)
    key <- matrix(
        c(rep(0, n_dna), rep(1, n_beads)), 
        nrow = 1, ncol = n_dna + n_beads, 
        dimnames = list("is_bead", ms[m]))
    key <- data.frame(key, check.names = FALSE)
    is_bead <- .get_bead_inds(x, key, assays[2])
    
    # subset specified assay
    cs <- assay(x, assays[1])
    es <- assay(x, assays[2])
    
    # get all events that should be removed later
    # including bead-bead and cell-bead doublets
    ths <- rowMins(es[bead_chs, is_bead])
    rmv <- colSums(es[bead_chs, ] > ths) == n_beads
    x$remove <- rmv
    
    # trim tails
    y <- es[bead_chs, is_bead]
    meds <- rowMedians(y)
    mads <- rowMads(y) * trim
    diff <- abs(y - meds)
    ex <- colAnys(diff > mads)
    is_bead[which(is_bead)[ex]] <- FALSE
    x$is_bead <- is_bead
    
    # get baselines (global mean)
    if (is.null(norm_to)) {
        bl <- rowMeans(cs[bead_chs, is_bead])
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
    smooth0 <- t(apply(cs[bead_chs, is_bead], 1, runmed, k, "constant"))
    
    # compute slopes (baseline versus smoothed bead intensitites)
    # & linearly interpolate slopes at non-bead events
    slopes <- colSums(smooth0 * bl) / colSums(smooth0 ^ 2)
    slopes <- approx(ts[is_bead], slopes, ts, rule = 2)$y
    
    # normalize raw bead intensities via multiplication with slopes
    cs <- sweep(cs, 2, slopes, "*")
    c <- ifelse(overwrite, assays[1], "normcounts")
    assay(x, c, FALSE) <- cs
    
    # smooth normalized beads
    y <- cs[bead_chs, is_bead]
    smooth <- t(apply(y, 1, runmed, k, "constant"))
    
    ps <- NULL
    if (plot) {
        ps <- list(
            scatter = .plot_bead_scatter(x, dna_chs, bead_chs, assays[2]), 
            lines = .plot_smooth_beads(smooth0, smooth, ts[is_bead]))
    }
    
    # do arcsinh-transformation on normalized counts
    if (transform) {
        e <- ifelse(overwrite, assays[2], "normexprs")
        x <- .transform(x, cofactor, ain = c, aout = e)
    }
    
    rownames(x) <- chs0
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
