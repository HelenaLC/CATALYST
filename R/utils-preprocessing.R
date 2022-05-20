#' @importFrom flowCore fr_append_cols flowSet
.match_ffs <- function(ffs, fix_chs) {
    chs <- lapply(ffs, colnames)
    if (!all(mapply(identical, chs[1], chs[-1]))) {
        warning(
            "Channel names don't match between FCS files; ", 
            "panel discrepancies will be fixed ", 
            "keeping ", fix_chs, " channels.")
        switch(
            fix_chs, 
            common = {
                chs_use <- Reduce(intersect, chs)
                mtx <- matrix(TRUE, length(chs_use), length(ffs))
                for (i in seq_along(ffs))
                    ffs[[i]] <- ffs[[i]][, chs_use]
            }, 
            all = {
                chs_use <- unique(unlist(chs))
                mtx <- vapply(chs, \(.) chs_use %in% ., logical(chs_use))
                for (i in seq_along(ffs)) {
                    chs2add <- setdiff(chs_use, chs[[i]])
                    if (length(chs2add) != 0) {
                        n <- nrow(ffs[[i]])
                        y <- vapply(chs2add, \(.) numeric(n), numeric(n))
                        ffs[[i]] <- fr_append_cols(ffs[[i]], y)
                    }
                    ffs[[i]] <- ffs[[i]][, chs_use]
                }
            }
        )
    } else mtx <- NULL
    fs <- flowSet(ffs)
    return(list(fs, mtx))
}
                      
#' @importFrom methods is
#' @importFrom flowCore isFCSfile read.FCS
.read_fs <- function(x, fix_chs, ...) {
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
        args <- list(...)
        defs <- c("transformation", "truncate_max_range")
        for (. in defs) if (is.null(args[[.]])) args[[.]] <- FALSE
        ffs <- lapply(fcs, \(fnm) do.call(read.FCS, c(fnm, args)))
        tmp <- .match_ffs(ffs, fix_chs)
        fs <- tmp[[1]]; mtx <- tmp[[2]]
    } else {
        stop("Invalid argument 'x'; should be either a flowSet",
            " or a character string specifying the path to", 
            " a directory containing a set of FCS files.")
    }
    if (!exists("mtx"))
        mtx <- NULL
    return(list(fs, mtx))
}

#' @importFrom SingleCellExperiment int_metadata<-
#' @importFrom SummarizedExperiment assay assay<- rowData
.transform <- function(x, cf, 
    ain = "counts", aout = "exprs",
    dir = c("forwards", "backwards")) {
    dir <- match.arg(dir)
    chs <- channels(x)
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
