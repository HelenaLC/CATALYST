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
#' a character string. If specified, outputs will be generated in this location. 
#' If NULL (the default), \code{normCytof} will return a \code{\link{flowFrame}}
#' of the normalized data (if \code{remove=FALSE}) or a \code{\link{flowSet}} 
#' containing normalized cells and beads (if \code{remove=TRUE}).
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
#' 
#' @return
#' if \code{out_path=NULL} (the default) a \code{\link{flowFrame}} of the 
#' normalized data (if \code{remove=FALSE}) or \code{\link{flowSet}} containing 
#' normalized cells and beads (if \code{remove=TRUE}). Else, a character of the 
#' location where output FCS files and plots have been generated.
#' 
#' @examples
#' path <- system.file("extdata", package="CATALYST")
#' fcsFiles <- list.files(path, "data_\\d+.fcs", full.names=TRUE)
#' raw_data <- concatFCS(fcsFiles)
#' normCytof(x = raw_data, y = "dvs")
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
    definition=function(x, y, out_path=NULL,
        remove_beads=TRUE, norm_to=NULL, k=500, trim=5) {
    
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
    message("Identifying beads...")
    key <- data.frame(matrix(c(0, 0, rep(1, n_beads)), ncol=2+n_beads,
        dimnames=list("is_bead", c(191, 193, bead_ms))), check.names=FALSE)
    re <- assignPrelim(x, key, verbose=FALSE)
    re <- estCutoffs(re, verbose=FALSE)
    re <- applyCutoffs(re)
    bead_inds <- bc_ids(re) == "is_bead"
    
    # get all events that should be removed later
    # including bead-bead and cell-bead doublets
    min_bead_ints <- matrixStats::colMins(es_t[bead_inds, bead_cols])
    remove <- apply(es_t[, bead_cols], 1, function(i) { 
        above_min <- sapply(1:n_beads, function(j) sum(i[j] > min_bead_ints[j]))
        sum(above_min) == n_beads })
    
    # trim tails
    meds <- matrixStats::colMedians(es_t[bead_inds, bead_cols])
    mads <- matrixStats::colMads(   es_t[bead_inds, bead_cols]) * trim
    ex <- matrixStats::colAnys(
        abs(t(es_t[bead_inds, bead_cols]) - meds) > mads)
    bead_inds[which(bead_inds)[ex]] <- FALSE

    # get slopes - baseline (global mean) versus smoothed bead intensitites
    # and linearly interpolate slopes at non-bead events
    message("Computing normalization factors...")
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
        re <- assignPrelim(norm_to, key, verbose=FALSE)
        re <- estCutoffs(re, verbose=FALSE)
        re <- applyCutoffs(re, sep_cutoffs=.3)
        bead_inds2 <- bc_ids(re) == "is_bead"
        
        meds <- matrixStats::colMedians(es2_t[bead_inds2, bead_cols])
        mads <- matrixStats::colMads(   es2_t[bead_inds2, bead_cols]) * trim
        ex <- matrixStats::colAnys(
            abs(t(es2_t[bead_inds2, bead_cols]) - meds) > mads)
        bead_inds2[which(bead_inds2)[ex]] <- FALSE
        
        bead_es <- es2[bead_inds2, bead_cols]
        bead_ts <- es2[bead_inds2, time_col]
    }
    baseline <- colMeans(bead_es)
    bead_slopes <- rowSums(bead_es*baseline) / rowSums(bead_es^2)
    slopes <- approx(bead_ts, bead_slopes, es[, time_col])$y
    
    # normalize and write FCS of normalized data
    normed_es <- cbind(
        es[,  c(time_col, lgth_col)], 
        es[, -c(time_col, lgth_col)]*slopes)
    outNormed(x, normed_es, remove_beads, remove, out_path)    

    # bead intensitites smoothed by conversion to local medians
    smoothed_beads <- data.frame(
        es[bead_inds, time_col], 
        sapply(bead_cols, function(i) 
            stats::runmed(es[bead_inds, i], k, "constant")))

    # normalize raw bead intensities via multiplication with slopes
    smoothed_normed_beads <- data.frame(
        sapply(c(time_col, bead_cols), function(i) 
            stats::runmed(normed_es[bead_inds, i], k, "constant")))
    colnames(smoothed_beads) <- colnames(smoothed_normed_beads) <- 
        chs[c(time_col, bead_cols)]
    
    p1 <- paste0(sprintf("%2.2f", sum(bead_inds)/nrow(es)*100), "%")
    p2 <- paste0(sprintf("%2.2f", sum(remove)   /nrow(es)*100), "%")
    t1 <- paste("Beads used for normalization (singlets only):", p1)
    t2 <- paste0("All events removed ",
        "(including bead-/cell-bead doublets): ", p2, "\n")
    ps <- c(
        plotBeads(es_t, bead_inds, bead_cols, dna_cols, TRUE, FALSE, TRUE),
        plotBeads(es_t, remove,    bead_cols, dna_cols, FALSE, TRUE, TRUE))
    
    if (!is.null(out_path)) {
        pdf(file.path(out_path, "beads_gated.pdf"), 
            width=n_beads*5, height=12.8)
        grid.arrange(grobs=ps, row=3, ncol=n_beads, 
            widths=rep(5, n_beads), heights=c(2.8, 5, 5),
            top=textGrob(paste(t1, t2, sep="\n"), just="right",
                gp=gpar(fontface="bold", fontsize=15)))
        dev.off()
    } else {
        grid.arrange(grobs=ps, row=3, ncol=n_beads, 
            widths=rep(5, n_beads), heights=c(3, 5, 5),
            top=textGrob(paste(t1, t2, sep="\n"), just="right",
                gp=gpar(fontface="bold", fontsize=15)))
    }
     
    if (!is.null(out_path)) {
        pdf(file.path(out_path, "beads_before_vs_after.pdf"), 
            width=15, height=12.5)
        grid.arrange(nrow=3, heights=c(6, .5, 6),
            plotSmoothed(smoothed_beads, "Smoothed beads"), 
            rectGrob(gp=gpar(fill="white", col="white")),
            plotSmoothed(smoothed_normed_beads, "Smoothed normalized beads"))
        dev.off()
    } else {
        grid.arrange(nrow=3, heights=c(6, .5, 6),
            plotSmoothed(smoothed_beads, "Smoothed beads"), 
            rectGrob(gp=grid::gpar(fill="white", col="white")),
            plotSmoothed(smoothed_normed_beads, "Smoothed normalized beads"))
    }
    })

# ------------------------------------------------------------------------------

#' @rdname normCytof
setMethod(f="normCytof",
    signature=signature(x="character"),
    definition=function(x, y, out_path=NULL,
        remove_beads=TRUE, norm_to=NULL, trim=5) {
        if (length(x) != 1) 
            stop("'x' should be a single character or flowFrame.")
        if (sum(flowCore::isFCSfile(x)) != 1) 
            stop(x, " is not a valid FCS file.")
        x <- flowCore::read.FCS(x)
        normCytof(x, y, out_path=NULL, remove_beads=TRUE, norm_to=NULL, trim=5)
    })
