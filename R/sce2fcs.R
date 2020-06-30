#' @rdname sce2fcs
#' @title SCE to \code{flowFrame/Set}
#' 
#' @description 
#' If \code{split_by = NULL}, the input SCE is converted to a 
#' \code{\link[flowCore:flowFrame-class]{flowFrame}}. Otherwise, 
#' it is split into a \code{\link[flowCore:flowSet-class]{flowSet}} 
#' by the specified \code{colData} column. 
#' Any cell metadata (\code{colData}) and dimension reductions 
#' available in the SCE may be dropped or propagated to the output.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param split_by NULL or a character string 
#'   specifying a \code{colData(x)} column to split by.
#' @param assay a character string specifying 
#'   which assay data to use; valid values are \code{assayNames(x)}. 
#'   When writing out FCS files, this should correspond to count-like data!
#' @param keep_cd,keep_dr logials specifying whether cell metadata 
#'   (stored in \code{colData(x)}) and dimension reductions 
#'   (stored in \code{reducedDims(x)}), respectively,
#'   should be kept or dropped.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @return 
#' a \code{\link[flowCore:flowFrame-class]{flowFrame}} 
#' if \code{split_by = NULL}; otherwise a 
#' \code{\link[flowCore:flowSet-class]{flowSet}}.
#' 
#' @examples 
#' # PREPROCESSING
#' data(sample_ff, sample_key)
#' sce <- prepData(sample_ff, by_time = FALSE)
#' sce <- assignPrelim(sce, sample_key, verbose = FALSE)
#' 
#' # split SCE by barcode population
#' fs <- sce2fcs(sce, split_by = "bc_id")
#' 
#' # do some spot checks
#' library(flowCore)
#' library(SingleCellExperiment)
#' 
#' length(fs) == nrow(sample_key)
#' all(fsApply(fs, nrow)[, 1] == table(sce$bc_id))
#' identical(t(exprs(fs[[1]])), assay(sce, "exprs")[, sce$bc_id == "A1"])
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
#' @importFrom flowCore flowFrame parameters
#' @importFrom matrixStats colMaxs
#' @importFrom methods is
#' @importFrom SingleCellExperiment int_colData int_metadata
#' @importFrom SummarizedExperiment assay assayNames rowData
#' @export 

sce2fcs <- function(x, split_by = NULL, 
    keep_cd = FALSE, keep_dr = FALSE, assay = "counts") {
    # check validity of input arguments
    args <- as.list(environment())
    .check_args_sce2fcs(args)
    
    # (optionally) split cells according to 'split_by'
    if (!is.null(split_by)) {
        cs <- split(seq_len(ncol(x)), factor(x[[split_by]]))
        l <- lapply(cs, function(i) x[, i])
    } else {
        cs <- list(seq_len(ncol(x)))
        l <- list(x)
    }
    
    # get numeric cell metadata columns to be added into assay data
    if (keep_cd) {
        cols_keep <- vapply(colData(x), function(u) 
            suppressWarnings(!all(is.na(as.numeric(as.character(u))))), 
            logical(1))
        for (i in which(cols_keep))
            x[[i]] <- as.numeric(x[[i]])
    } else cols_keep <- FALSE
    
    # construct 'flowFrame' identifiers
    if (!is.null(split_by)) {
        fix_name <- vapply(seq_along(l), function(i) 
            !is.na(suppressWarnings(as.numeric(names(l)[i]))),
            logical(1))
        prefix <- ifelse(fix_name, paste0(split_by, "_"), "")
        ids <- paste0(prefix, names(l))
    } else ids <- NULL
    
    # get internal cell metadata (optionally) including reduced dimensions
    icd <- as.matrix(data.frame(int_colData(x)))
    is_dr <- grepl("reducedDims", colnames(icd))
    pat <- "reducedDims.([[:alpha:]]+).([0-9]+)"
    colnames(icd) <- gsub(pat, "\\1\\2", colnames(icd))
    if (!keep_dr) icd <- icd[, !is_dr]
    
    # construct expression matrix including all cell metadata
    y <- t(assay(x, assay))
    cd <- colData(x)[cols_keep]
    y <- cbind(icd, y, as.matrix(cd))
    
    is <- seq_along(l)
    names(is) <- ids
    ffs <- lapply(is, function(i) {
        # construct preliminary 'flowFrame'
        z <- y[cs[[i]], , drop = FALSE]
        ff <- flowFrame(z)
        # update 'AnnotatedDataFrame' of parameters
        ps <- parameters(ff)
        rs <- 2^(ceiling(log2(colMaxs(z))))
        m <- match(rownames(x), colnames(z), nomatch = 0)
        ps$name[m] <- channels(x)
        ps$desc[-m] <- NA
        ps$range <- rs
        ps$minRange <- 0
        ps$maxRange <- rs-1
        colnames(z) <- ps$name
        # update description list
        ds <- int_metadata(x)$description
        if (is.null(ds)) 
            ds <- list()
        if (is.null(ds$`$CYT`))
            ds$`$CYT` <- "cytof2"
        ds$GUID <- ids[i]
        ds$transformation <- "custom"
        ds$`$DATE` <- format(Sys.Date(), "%d-%b-%Y")
        for (i in seq_len(nrow(ps))) {
            ds[[sprintf("$P%sN", i)]] <- ps$name[i]
            ds[[sprintf("$P%sS", i)]] <- ps$desc[i]
            ds[[sprintf("$P%sR", i)]] <- paste(rs[i])
            ds[[sprintf("flowCore_$P%sRmin", i)]] <- "0"
            ds[[sprintf("flowCore_$P%sRmax", i)]] <- paste(rs[i]-1)
        }
        # construct final 'flowFrame'
        attr(z, "ranges") <- rs-1
        flowFrame(z, ps, ds)
    })
    if (length(ffs) == 1) ffs[[1]] else flowSet(ffs)
}
