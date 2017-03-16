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
#' @param x a character string of the FCS file to be normalized.
#' @param y \code{"dvs"} or \code{"beta"} or a numeric vector of bead masses.
#' @param gate \code{"auto"} or \code{"manual"}. 
#' Specifies how beads should be gated.
#' @param out_path 
#' a character string. If specified, outputs will be generated in this location. 
#' If NULL (the default), \code{normCytof} will return a \code{\link{flowFrame}}
#' of the normalized data (if \code{remove=FALSE}) or a \code{\link{flowSet}} 
#' containing normalized cells and beads (if \code{remove=TRUE}).
#' @param remove_beads logical. If TRUE (the default) normalized cells and beads 
#' will be returned separately.
#' @param norm_to a character string of an FCS file from which baseline values
#' should be computed and to which the input data should be normalized.
#' 
#' @return
#' if \code{out_path=NULL} (the default) a \code{\link{flowFrame}} of the 
#' normalized data (if \code{remove=FALSE}) or \code{\link{flowSet}} containing 
#' normalized cells and beads (if \code{remove=TRUE}). Else, a character of the 
#' location where output FCS files and plots have been generated.
#' 
#' @examples
#'
#' @references 
#' Finck, R. et al. (2013).
#' Normalization of mass cytometry data with bead standards.
#' \emph{Cytometry A} \bold{83A}, 483–494.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2 grid gridExtra
#' @importFrom RColorBrewer brewer.pal

# ------------------------------------------------------------------------------

x <- "/Users/HLC/Dropbox/spillover/concatenation test/20150918_MDMenz_concat.fcs"
normCytof(x, "dvs", out_path="/Users/HLC/Dropbox/spillover/concatenation test/")

setMethod(f="normCytof", 
    signature=signature(x="character"), 
    definition=function(x, y, gate="auto", 
        out_path=NULL, remove_beads=TRUE, norm_to=NULL) {
    
    ff <- flowCore::read.FCS(x)
    es <- flowCore::exprs(ff)
    es_t <- asinh(es/5)
    chs <- flowCore::colnames(ff)
    ms <- gsub("[[:alpha:][:punct:]]", "", chs)
    ms <- ms[!is.na(as.numeric(ms))]
    
    # find time, length, DNA and bead channels
    time_col <- grep("time", chs, TRUE)
    length_col <- grep("length", chs, TRUE)
    dna_cols <- grep("Ir191|Ir193", chs, TRUE)
    tmp <- get_bead_cols(chs, y)
    bead_cols <- tmp[[1]]
    bead_ms <- tmp[[2]]
    n_beads <- length(bead_ms)

    # identify bead singlets
    message("Identifying beads...")
    key <- data.frame(matrix(c(0, 0, rep(1, n_beads)), ncol=2+n_beads,
        dimnames=list("is_bead", c(191, 193, bead_ms))), check.names=FALSE)
    re <- assignPrelim(ff, key, verbose=FALSE)
    re <- estCutoffs(re, verbose=FALSE)
    re <- applyCutoffs(re, sep_cutoffs=0.5)
    bead_inds <- bc_ids(re) == "is_bead"
    ex <- sapply(bead_cols, function(k) 
        abs(es_t[bead_inds, k] - median(es_t[bead_inds, k])) > 
            2*mad(es_t[bead_inds, k]))
    ex <- apply(ex, 1, any)
    bead_inds[ex] <- FALSE
    
    # get all events that should be removed later
    # including bead-bead and cell-bead doublets
    min_bead_ints <- apply(es_t[bead_inds, bead_cols], 2, min)
    remove <- apply(es_t[, bead_cols], 1, function(k) { 
        above_min <- sapply(1:n_beads, function(i) sum(k[i] > min_bead_ints[i]))
        sum(above_min) == n_beads })
    
    # get slopes - baseline (global mean) versus smoothed bead intensitites
    if (is.null(norm_to)) {
        message("Computing normalization factors...")
        bead_es <- es[bead_inds, bead_cols]
    } else {
        bead_ff <- flowCore::read.FCS(norm_to)
        bead_es <- flowCore::exprs(bead_ff[, bead_cols])
    }
    baseline <- apply(bead_es, 2, mean)
    bead_slopes <- apply(bead_es*baseline, 1, sum) / apply(bead_es^2, 1, sum)
    
    # linearly interpolate slopes
    times <- es[, time_col]
    slopes <- approx(times[bead_inds], bead_slopes, times)$y
    
    # normalize 
    normed_es <- cbind(
        es[,  c(time_col, length_col)], 
        es[, -c(time_col, length_col)]*slopes)

    if (remove_beads) {
        cells <- new("flowFrame",
            exprs=normed_es[!remove, ],
            parameters=flowCore::parameters(ff),
            description=flowCore::description(ff))
        beads <- new("flowFrame",
            exprs=normed_es[remove, ],
            parameters=flowCore::parameters(ff),
            description=flowCore::description(ff))
        if (is.null(out_path)) {
            flowSet(cells, beads)
        } else {
            out_nm <-  file.path(gsub(".fcs", "", x))
            suppressWarnings(flowCore::write.FCS(cells, 
                paste0(out_nm, "_normalized.fcs")))
            suppressWarnings(flowCore::write.FCS(beads, 
                paste0(out_nm, "_beads.fcs")))
        }
    } else {
        normed <- new("flowFrame",
            exprs=normed_es,
            parameters=flowCore::parameters(ff),
            description=flowCore::description(ff))
        if (is.null(out_path)) {
            normed
        } else {
            suppressWarnings(flowCore::write.FCS(normed, 
                paste0(file.path(gsub(".fcs", "", x)), "_normalized.fcs")))
        }
    }

    # bead intensitites smoothed by conversion to local medians
    smoothed_beads <- data.frame(
        es[bead_inds, time_col], 
        sapply(bead_cols, function(k) 
            runmed(es[bead_inds, k], 501, "constant")))

    # normalize raw bead intensities via multiplication with slopes
    smoothed_normed_beads <- data.frame(
        sapply(c(time_col, bead_cols), function(k) 
            runmed(normed_es[bead_inds, k], 501, "constant")))
    colnames(smoothed_beads) <- colnames(smoothed_normed_beads) <- 
        chs[c(time_col, bead_cols)]
    
    m <- matrix(1:(2*n_beads), nrow=2)
    m <- m[nrow(m):1, ]
    pdf(file.path(out_path, "beads.pdf"), width=n_beads*3)
    gridExtra::grid.arrange(nrow=3, heights=c(2.5,4,4), 
        layout_matrix=rbind(m, (2*n_beads+1):(3*n_beads)), grobs=c(
            plotBeads(es_t, bead_inds, bead_cols, dna_cols, TRUE, FALSE, TRUE),
            plotBeads(es_t, remove,    bead_cols, dna_cols, FALSE, TRUE, TRUE)))
    dev.off()
   
    pdf(file.path(out_path, "beads_before_vs_after.pdf"), width=15, height=12.5)
    grid.arrange(nrow=3, heights=c(6,.5,6),
        plotSmoothedBeads(smoothed_beads, "Smoothed beads"), 
        rectGrob(gp=gpar(fill="white", col="white")),
        plotSmoothedBeads(smoothed_normed_beads, "Smoothed normalized beads"))
    dev.off()
    })

