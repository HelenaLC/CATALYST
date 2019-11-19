#' @importFrom Biobase pData
#' @importFrom flowCore isFCSfile read.flowSet fsApply exprs parameters description
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export

fcs2sce <- function(x) {
    x <- c(
        "~/Dropbox/spillover new/bead based compensation/bead replicates/160616_beads.fcs",
        "~/Dropbox/spillover new/bead based compensation/bead replicates/160805_beads.fcs")
    stopifnot(isFCSfile(x))
    x <- read.flowSet(x)
    y <- fsApply(x, exprs)
    ps <- fsApply(x, parameters)
    ds <- fsApply(x, description)
    
    x <- c("~/Dropbox/spillover new/bead based compensation/bead replicates/160616_beads.fcs")
    x <- read.FCS(x)
    es <- exprs(x)
    chs <- colnames(es)
    ms <- .get_ms_from_chs(chs)
    mets <- .get_mets_from_chs(chs)
    
    cd <- NULL
    rd <- data.frame(pData(parameters(x)), row.names = NULL)
    md <- description(x)
    es <- t(es)
    rownames(es) <- colnames(es) <- NULL
    sce <- SingleCellExperiment(
        assays = list(counts = es),
        rowData = rd, metadata = md)
    
    sub <- sce[, sample(ncol(sce), 4e3)]
    assay(sub, "logcounts") <- asinh(assay(sub)/5)
    sub <- runUMAP(sub)
    sub <- runTSNE(sub)
    plotUMAP(sub)
    plotTSNE(sub)
}
