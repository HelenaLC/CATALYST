#' @rdname fcs2sce
#' @title `SingleCellExperiment` from FCS file(s)
#' 
#' @description Concatenates input FCS files to a `SingleCellExperiment`.
#' 
#' @param x character vector of FCS file names.
#' @param by_time logical. Should files be ordered by acquisition time?
#' @param file_no logical. Should a file number parameter be added?
#' 
#' @author Helena L. Crowell
#' 
#' @importFrom Biobase pData
#' @importFrom flowCore isFCSfile read.flowSet fsApply exprs parameters description keyword
#' @importFrom matrixStats colMaxs
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export

fcs2sce <- function(x, by_time = TRUE, file_no = FALSE) {
    stopifnot(
        is(x, "flowFrame") || is(x, "flowSet") 
            || is(x, "character") && isFCSfile(x),
        is.logical(by_time), length(by_time) == 1,
        is.logical(file_no), length(file_no) == 1)
    
    if (is(x, "character")) {
        fns <- x
        fs <- read.flowSet(x, 
            transformation = FALSE, 
            truncate_max_range = FALSE)
    } else {
        if (is(x, "flowFrame"))
            x <- flowSet(x)
        fns <- fsApply(x, identifier)
        fs <- x
    }
    
    if (by_time) {
        ts <- keyword(fs, "$BTIM")
        if(any(vapply(ts, is.null, logical(1)))) {
            message("Not all samples contain information on their",
                " acquisition time; ignoring argument 'by_time'.",
                " Samples will be kept in their original order.")
        } else {
            fs <- fs[order(ts)]
        }
    }
    
    y <- fsApply(fs, exprs)
    ns <- fsApply(fs, nrow)
    ps <- fsApply(fs, parameters)
    ds <- fsApply(fs, description)
    
    if (file_no) y[, "FileNo"] <- rep.int(seq_along(fs), ns)
    
    t <- grep("time", colnames(fs), ignore.case = TRUE)
    t0 <- c(1, cumsum(ns) + 1)
    tx <- t0[-1] - 1
    for (i in seq_along(fs)[-1]) {
        idx <- seq(t0[i], tx[i])
        y[idx, t] <- y[idx, t] + y[tx[i - 1], t]
    }
    
    # construct 'SingleCellExperiment'
    cd <- data.frame(file_name = rep.int(basename(fns), ns))
    rd <- data.frame(pData(ps[[1]]), row.names = 1)
    rd$maxRange <- ceiling(colMaxs(fsApply(fs, function(u) colMaxs(exprs(u)))))
    rd$range <- rd$maxRange + 1
    
    SingleCellExperiment(
        assays = list(counts = t(y)), 
        colData = cd, rowData = rd, 
        metadata = list(description = ds))
}
