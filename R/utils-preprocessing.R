
#' @importFrom methods is
#' @importFrom flowCore isFCSfile read.flowSet
.read_fs <- function(x) {
    if (is(x, "flowFrame")) {
        fs <- flowSet(x)
    } else if (is.list(x)) {
        stopifnot(all(vapply(x, function(u) 
            is(u, "flowFrame"), logical(1))))
        fs <- flowSet(x)
    } else if (is(x, "flowSet")) {
        fs <- x
    } else if (is.character(x)) {
        if (length(x) == 1 && dir.exists(x)) {
            fcs <- list.files(x, 
                pattern = ".fcs$", 
                full.names = TRUE,
                ignore.case = TRUE)
            if (length(fcs) == 0) 
                stop("Couldn't find any FCS files",
                    " in the specified directory.")
        } else {
            fcs <- x  
        }
        stopifnot(all(vapply(fcs, isFCSfile, logical(1))))
        fs <- read.flowSet(fcs,
            transformation = FALSE, 
            truncate_max_range = FALSE)
    } else {
        stop("Invalid argument 'x'; should be either a flowSet",
            " or a character string specifying the path to", 
            " a directory containing a set of FCS files.")
    }
    return(fs)
}

#' @importFrom SingleCellExperiment int_metadata<-
#' @importFrom SummarizedExperiment assay assay<- rowData
.transform <- function(x, cf, 
    ain = "counts", aout = "exprs",
    dir = c("forwards", "backwards")) {
    dir <- match.arg(dir)
    chs <- rowData(x)$channel_name
    stopifnot(is.numeric(cf), cf > 0)
    if (length(cf) == 1) {
        int_metadata(x)$cofactor <- cf
        cf <- rep(cf, nrow(x))
    } else {
        stopifnot(!is.null(names(cf)), chs %in% names(cf))
        cf <- cf[match(chs, names(cf))]
        int_metadata(x)$cofactor <- cf
    }
    if (dir == "forwards") {
        fun <- asinh
        op <- "/"
    } else {
        fun <- sinh
        op <- "*"
    }
    y <- assay(x, ain)
    y <- fun(sweep(y, 1, cf, op))
    assay(x, aout, FALSE) <- y
    return(x)
}
