#' @rdname sce2fcs
#' @title SCE to \code{flowFrame/Set}
#' 
#' @description If \code{split_by = NULL}, the input SCE is converted to a 
#' \code{\link[flowCore]{flowFrame}}. Otherwise, the SCE is split into a 
#' \code{\link[flowCore]{flowSet}} by the specified \code{colData} column.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param split_by NULL or a character string 
#'   specifying a \code{colData(x)} column to split by.
#' @param assay a character string specifying which assay data to use.
#'   Valid values are \code{assayNames(x)}.
#' @param keep_cd,keep_dr logials specifying whether cell metadata 
#'   (stored in \code{colData(x)}) and dimension reductions 
#'   (stored in \code{reducedDims(x)}), respectively,
#'   should be kept or dropped.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @return a \code{\link[flowCore]{flowFrame}} if \code{split_by = NULL};
#' a \code{\link[flowCore]{flowSet}}, if \code{split_by} corresponds to a
#' \code{colData} column of \code{x}.
#' 
#' @examples 
#' # PREPROCESSING
#' data(sample_ff, sample_key)
#' sce <- fcs2sce(sample_ff, by_time = FALSE)
#' sce <- assignPrelim(sce, sample_key, verbose = FALSE)
#' 
#' # split SCE by barcode population
#' fs <- sce2fcs(sce, split_by = "bc_id", assay = "exprs")
#' 
#' # do some spot checks
#' library(flowCore)
#' length(fs) == nrow(sample_key)
#' all(fsApply(fs, nrow)[, 1] == table(sce$bc_id))
#' identical(t(exprs(fs$A1)), assay(sce, "exprs")[, sce$bc_id == "A1"])
#' 
#' # DIFFERENTIAL ANALYSIS
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' sce <- cluster(sce, verbose = FALSE)
#' 
#' # split by 20 metacluster populations
#' sce$meta20 <- cluster_ids(sce, "meta20")
#' fs <- sce2fcs(sce, split_by = "meta20", assay = "exprs")
#' all(fsApply(fs, nrow)[, 1] == table(sce$meta20))
#' 
#' @importFrom flowCore flowFrame
#' @importFrom methods as is
#' @importFrom SummarizedExperiment assay assayNames
#' @export 

sce2fcs <- function(x, 
    split_by = NULL, assay = "exprs",
    keep_cd = FALSE, keep_dr = FALSE) {
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"), is.null(split_by) 
        || is.character(split_by) & !is.null(x[[split_by]]),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        is.logical(keep_cd), length(keep_cd) == 1,
        is.logical(keep_dr), length(keep_dr) == 1)
    if (keep_dr && length(drs <- reducedDims(x)) != 0) {
        # concatenate dimension reductions
        n_dims <- vapply(drs, ncol, numeric(1))
        drs <- do.call("cbind", lapply(drs, data.frame))
        colnames(drs) <- paste(sep = "_",
            rep.int(reducedDimNames(x), n_dims),
            unlist(lapply(n_dims, seq_len)))
    } else keep_dr <- FALSE
    if (!is.null(split_by)) {
        cs <- split(seq_len(ncol(x)), factor(x[[split_by]]))
        l <- lapply(cs, function(i) x[, i])
    } else {
        cs <- list(seq_len(ncol(x)))
        l <- list(x)
    }
    ds <- metadata(x)
    ffs <- lapply(seq_along(l), function(i) {
        y <- t(assay(l[[i]], assay))
        if (keep_cd) {
            cols_keep <- vapply(colData(x), function(u) 
                suppressWarnings(!all(is.na(as.numeric(as.character(u))))), 
                logical(1))
            cd <- as.matrix(colData(x)[cs[[i]], cols_keep, drop = FALSE])
            y <- cbind(y, cd)
        }
        if (keep_dr) y <- cbind(y, as.matrix(drs[cs[[i]], ]))
        ff <- flowFrame(y)
        ps <- parameters(ff)
        for (i in seq_len(nrow(ps))) {
            ds[[sprintf("$P%sN", i)]] <- as.character(ps$name[i])
            ds[[sprintf("$P%sS", i)]] <- as.character(ps$desc[i])
        }
        description(ff) <- ds
        return(ff)
    })
    fs <- as(ffs, "flowSet")
    for (i in seq_along(fs)) {
        id <- paste(split_by, names(l)[i], sep = ".")
        id <- ifelse(length(id) == 0, NA, id)
        description(fs[[i]])[c("GUID", "ORIGINALGUID")] <- 
            identifier(fs[[i]]) <- id
    }
    if (length(fs) == 1) fs[[1]] else fs
}
