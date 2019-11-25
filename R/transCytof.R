#' @rdname transform
#' @title ...
#' 
#' @description ...
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param method character or function.
#' @param cf numeric cofactor for arcsinh-transformation.
#' @param assay character. Valid values are \code{assayNames}.
#' @param id character. ID to be given to newly generated assay slot.
#' @param overwrite logical. Should assay slot \code{id} be overwritten?
#' @param ... optional arguments passed to \code{method}.
#' 
#' @examples 
#' data(raw_data)
#' sce <- fcs2sce(raw_data)
#' 
#' @author Helena L. Crowell
#' 
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay assay<- assayNames
#' @export

transform <- function(x, method = "arcsinh", cf = 10, 
    assay = "counts", id = "exprs", overwrite = FALSE, ...) {
    stopifnot(is(x, "SingleCellExperiment"),
        is.numeric(cf), length(cf) == 1, cf > 0,
        is.character(assay), length(assay) == 1, assay %in% assayNames(x))
    if (id %in% assayNames(x) & !overwrite)
        stop("Assay slot ", dQuote(id), " already exists;",
            " specify another 'id' or use 'overwrite = TRUE'.")
    if (is.function(method)) {
        t <- method
    } else {
        method <- match.arg(method)
        t <- switch(method, "arcsinh" = function(u, ...) asinh(u/cf))
    }
    ms <- .get_ms_from_chs(rownames(x))
    y <- assay(x, assay)
    y[!is.na(ms), ] <- t(y[!is.na(ms), ], ...)
    assay(x, id) <- y
    return(x)
}
