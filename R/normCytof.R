# ==============================================================================
# normalization with bead standards
# ------------------------------------------------------------------------------

#' @rdname normCytof
#' @title Bead-based normalization
#' 
#' @description 
#' an implementation of Finck et al.'s normalization of mass cytometry data 
#' using bead standards with automated bead gating.
#'
#' @param x 
#' a \code{\link{flowFrame}} or character of the FCS file to be normalized.
#' @param y 
#' \code{"dvs"} (for bead masses 140, 151, 153, 165, 175) or \code{"beta"} 
#' (for bead masses 139, 141, 159, 169, 175) or a numeric vector of bead masses.
#' @param out_path 
#' a character string. If specified, outputs will be generated 
#' in this location. If NULL (the default), \code{normCytof} will return 
#' a \code{\link{flowFrame}} of the normalized data (if \code{remove=FALSE}) 
#' or a \code{\link{flowSet}} containing normalized cells and beads 
#' (if \code{remove=TRUE}).
#' @param remove_beads 
#' logical. If TRUE (the default) beads will be removed and normalized cells and
#' beads returned separately.
#' @param norm_to 
#' a \code{\link{flowFrame}} or character of an FCS file from which baseline 
#' values should be computed and to which the input data should be normalized.
#' @param k integer width of the median window used for bead smoothing.
#' @param trim 
#' a single non-negative numeric. A \emph{median} +/- ... \emph{mad} rule is
#' applied to the preliminary population of bead events to remove bead-bead 
#' doublets and low signal beads prior to estimating normalization factors.
#' @param verbose  
#' logical. Should extra information on progress be reported? Defaults to TRUE.
#' 
#' @return
#' if \code{out_path=NULL} (the default) a \code{\link{flowFrame}} of the 
#' normalized data (if \code{remove=FALSE}) or \code{\link{flowSet}} containing 
#' normalized cells and beads (if \code{remove=TRUE}). Else, a character of the 
#' location where output FCS files and plots have been generated.
#' 
#' @examples
#' data(raw_data)
#' ff <- concatFCS(raw_data)
#' normCytof(x = ff, y = "dvs", k=300)
#'
#' @references 
#' Finck, R. et al. (2013).
#' Normalization of mass cytometry data with bead standards.
#' \emph{Cytometry A} \bold{83A}, 483-494.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @import matrixStats
#' @importFrom flowCore colnames exprs flowFrame flowSet read.FCS write.FCS
#' @importFrom grDevices pdf dev.off
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats approx mad runmed

# ------------------------------------------------------------------------------

setMethod(f="normCytof", 
    signature=signature(x="flowFrame"), 
    definition=function(x, y, out_path=NULL, remove_beads=TRUE, norm_to=NULL, 
        k=500, trim=5, verbose=TRUE) {
    
    # assure width of median window is odd
    if (k %% 2 == 0) k <- k + 1
    
    es <- flowCore::exprs(x)
    es_t <- asinh(es/5)
    chs <- flowCore::colnames(x)
    ms <- gsub("[[:alpha:][:punct:]]", "", chs)
    
    # find time, length, DNA and bead channels
    time_col <- grep("time",        chs, ignore.case=TRUE)
    lgth_col <- grep("length",      chs, ignore.case=TRUE)
    dna_cols <- grep("Ir191|Ir193", chs, ignore.case=TRUE)
    bead_cols <- get_bead_cols(chs, y)
    bead_ms <- ms[bead_cols]
    n_beads <- length(bead_ms)
    ms <- ms[!is.na(as.numeric(ms))]

    # identify bead singlets
    if (verbose) message("Identifying beads...")
    key <- data.frame(matrix(c(0, 0, rep(1, n_beads)), ncol=2+n_beads,
        dimnames=list("is_bead", c(191, 193, bead_ms))), check.names=FALSE)
    bead_inds <- get_bead_inds(x, key)
    
    # get all events that should be removed later
    # including bead-bead and cell-bead doublets
    min_bead_ints <- matrixStats::colMins(es_t[bead_inds, bead_cols])
    remove <- apply(es_t[, bead_cols], 1, function(i) { 
        above_min <- vapply(seq_len(n_beads), 
            function(j) sum(i[j] > min_bead_ints[j]), numeric(1))
        sum(above_min) == n_beads 
    })
    
    # trim tails
    bead_inds <- update_bead_inds(es_t, bead_inds, bead_cols, trim)

    # get slopes - baseline (global mean) versus smoothed bead intensitites
    # and linearly interpolate slopes at non-bead events
    if (verbose) message("Computing normalization factors...")
    if (is.null(norm_to)) {
        bead_es <- es[bead_inds, bead_cols]
        bead_ts <- es[bead_inds, time_col]
    } else {
        if (is.character(norm_to)) {
            if (length(norm_to) != 1) 
                stop("'norm_to' should be a single character or flowFrame.")
            if (sum(flowCore::isFCSfile(norm_to)) != 1) 
                stop(norm_to, " is not a valid FCS file.")
            norm_to <- flowCore::read.FCS(norm_to)
        }
        es2 <- flowCore::exprs(norm_to)
        es2_t <- asinh(es2/5)
        bead_inds2 <- get_bead_inds(norm_to, key)
        bead_inds2 <- update_bead_inds(es2_t, bead_inds2, bead_cols, trim)
        
        bead_es <- es2[bead_inds2, bead_cols]
        bead_ts <- es2[bead_inds2, time_col]
    }
    baseline <- colMeans(bead_es)
    bead_slopes <- rowSums(bead_es*baseline) / rowSums(bead_es^2)
    slopes <- approx(bead_ts, bead_slopes, es[, time_col])$y
    
    # normalize data
    normed_es <- cbind(
        es[,  c(time_col, lgth_col)], 
        es[, -c(time_col, lgth_col)]*slopes)

    # bead intensitites smoothed by conversion to local medians
    smoothed_beads <- data.frame(
        es[bead_inds, time_col], 
        vapply(bead_cols, function(i) 
            stats::runmed(es[bead_inds, i], k, "constant"),
            numeric(sum(bead_inds))))

    # normalize raw bead intensities via multiplication with slopes
    smoothed_normed_beads <- data.frame(
        es[bead_inds, time_col], 
        vapply(bead_cols, function(i) 
            stats::runmed(normed_es[bead_inds, i], k, "constant"),
            numeric(sum(bead_inds))))
    colnames(smoothed_beads) <- colnames(smoothed_normed_beads) <- 
        chs[c(time_col, bead_cols)]
    
    outPlots(es_t, bead_inds, remove, bead_cols, dna_cols,
        smoothed_beads, smoothed_normed_beads, out_path)
    outNormed(x, normed_es, remove_beads, remove, out_path)
    })

# ------------------------------------------------------------------------------

#' @rdname normCytof
setMethod(f="normCytof",
    signature=signature(x="character"),
    definition=function(x, y, out_path=NULL, remove_beads=TRUE, norm_to=NULL, 
        k=500, trim=5, verbose=TRUE) {
        if (length(x) != 1) 
            stop("'x' should be a single character or flowFrame.")
        if (sum(flowCore::isFCSfile(x)) != 1) 
            stop(x, " is not a valid FCS file.")
        x <- flowCore::read.FCS(x)
        normCytof(x, y, out_path=NULL, remove_beads=TRUE, norm_to=NULL, 
            k=500, trim=5, verbose=TRUE)
    })
