#' @rdname fcs2sce
#' @title FCS file & \code{flowFrame/-Set} 
#'   conversion to \code{SingleCellExperiment}
#' 
#' @description Converts (and concatenates, if multiple samples are supplied) 
#'   input FCS file(s), \code{flowFrame}(s), or a \code{flowSet} to an object 
#'   of the \code{SingleCellExperiment} class.
#' 
#' @param x character vector of FCS file names, 
#'   \code{\link[flowCore]{flowFrame}}, or \code{\link[flowCore]{flowSet}}.
#' @param by_time logical. Should files be ordered 
#'   by acquisition time? (see details)
#' @param file_no logical. Should a file number parameter be stored
#'   in the \code{int_colData} of the output SCE?
#' @param transform logical. Specifies whether an arcsinh-transformation
#'   with cofactor \code{cofactor} should be performed, 
#'   in which case expression values (transformed counts) 
#'   will be stored in \code{assay(x, "exprs")}.
#' @param cofactor numeric cofactor to use for optional 
#'   arcsinh-transformation if \code{transform = TRUE}.
#' 
#' @details By default, non-mass channels (e.g., time, event lengths)
#'   will be removed from the output SCE's assay data and instead
#'   stored in the object's internal cell metadata (\code{int_colData})
#'   to assure these data are not subject to transformations
#'   or other computations applied to the assay data.
#' 
#'   When \code{x} contains more than one sample, \code{fcs2sce} will 
#'   concatenate cells into a single \code{SingleCellExperiment} object.
#'   Note that cells will hereby be order by \code{"Time"}, regardless 
#'   of whether \code{by_time = TRUE} or \code{FALSE}. 
#'   Instead, \code{by_time} determines the sample (not cell!) order; 
#'   i.e., whether samples should be kept in their original order,
#'   or should be re-ordered according to their acquision time 
#'   stored in \code{keyword(flowSet, "$BTIM")}.
#'   
#' @return a \code{SingleCellExperiment}.
#' 
#' @examples
#' data(raw_data)
#' (sce <- fcs2sce(raw_data))
#' 
#' library(SingleCellExperiment)
#' assayNames(sce)         # view available assays
#' names(int_colData(sce)) # non-mass channels
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @importFrom Biobase pData
#' @importFrom flowCore isFCSfile read.flowSet 
#'   fsApply exprs parameters description keyword
#' @importFrom matrixStats colMaxs
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment 
#'   SingleCellExperiment int_colData<- int_metadata<-
#' @importFrom SummarizedExperiment assay assay<-
#' @export

fcs2sce <- function(x, 
    by_time = TRUE, file_no = FALSE, 
    transform = TRUE, cofactor = 5) {
    stopifnot(
        is(x, "flowFrame") || is(x, "flowSet") 
        || is(x, "character") && isFCSfile(x),
        is.logical(by_time), length(by_time) == 1,
        is.logical(file_no), length(file_no) == 1,
        is.logical(transform), length(transform) == 1,
        is.numeric(cofactor), length(cofactor) == 1, cofactor > 0)
    
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
    
    if (by_time & length(fs) > 1) {
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
    
    # fix event times
    t <- grep("time", colnames(fs), ignore.case = TRUE)
    if (length(t) != 0) {
        t0 <- c(1, cumsum(ns) + 1)
        tx <- t0[-1] - 1
        for (i in seq_along(fs)[-1]) {
            idx <- seq(t0[i], tx[i])
            y[idx, t] <- y[idx, t] + y[tx[i - 1], t]
        }
    }
    
    # construct SCE
    cd <- data.frame(file_name = rep.int(basename(fns), ns))
    rd <- data.frame(pData(ps[[1]]), row.names = 1)
    #rd$maxRange <- ceiling(colMaxs(fsApply(fs, function(u) colMaxs(exprs(u)))))
    #rd$range <- rd$maxRange + 1
    z <- t(y); rownames(z) <- unname(colnames(y))
    sce <- SingleCellExperiment(
        assays = list(counts = z),
        colData = cd, rowData = rd, 
        metadata = list(description = ds))
    int_metadata(sce)$description <- ds
    
    # grep non-mass channels, store in internal 
    # event metadata, and exclude from count matrix
    is_mass <- !is.na(.get_ms_from_chs(colnames(y)))
    foo <- DataFrame(matrix(vector(), nrow = nrow(y)))
    int_cd <- DataFrame(y[, !is_mass], check.names = FALSE)
    colnames(int_cd) <- colnames(y)[!is_mass]
    int_cd$reducedDims <- int_cd$altExps <- foo
    # (optionally) add file number & return SCE
    if (file_no) int_cd$file_no <- rep.int(seq_along(fs), ns)
    int_colData(sce) <- int_cd
    sce <- sce[is_mass, ]
    
    # (optionally) apply arcsinh-transformation
    if (transform) {
        metadata(sce)$cofactor <- cofactor
        assay(sce, "exprs", FALSE) <- asinh(assay(sce) / cofactor)
    }
    return(sce)
}
