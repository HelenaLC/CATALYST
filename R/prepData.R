#' @rdname prepData
#' @title Data preparation
#' 
#' @param x a \code{flowSet} holding all samples 
#'   or a path to a set of FCS files.
#' @param panel a data.frame containing, for each channel, 
#'   its column name in the input data, targeted protein marker,
#'   and (optionally) class ("type", "state", or "none").
#'   If `panel` is unspecified, it will be constructed 
#'   from the first input sample via \code{\link{guessPanel}}.
#' @param md a table with column describing the experiment.
#'   An exemplary metadata table could look as follows: \itemize{
#'   \item\code{file_name}: the FCS file name
#'   \item\code{sample_id}: a unique sample identifier
#'   \item\code{patient_id}: the patient ID
#'   \item\code{condition}: brief sample description 
#'     (e.g. reference/stimulated, healthy/diseased)}
#'   If `md` is unspecified, the \code{flowFrame/Set} 
#'   \code{\link[flowCore]{identifier}}(s) will be used
#'   as sample IDs with no additional metadata factors.
#' @param features a logical vector, numeric vector of column indices,
#'   or character vector of channel names. Specified which column to keep 
#'   from the input data. Defaults to the channels listed in the input panel.
#' @param transform logical. Specifies whether an arcsinh-transformation with 
#'   cofactor cofactor should be performed, in which case expression values 
#'   (transformed counts) will be stored in assay(x, "exprs").
#' @param cofactor numeric cofactor(s) to use for optional 
#'   arcsinh-transformation when \code{transform = TRUE};
#'   single value or a vector with channels as names.
#' @param panel_cols a names list specifying the column names of \code{panel}
#'   that contain the channel names, targeted protein markers, and (optionally) 
#'   marker classes. When only some `panel_cols` deviate from the defaults,
#'   specifying only these is sufficient.
#' @param md_cols a named list specifying the column names of \code{md}
#'   that contain the FCS file names, sample IDs, and factors of interest
#'   (batch, condition, treatment etc.). When only some `md_cols` deviate 
#'   from the defaults, specifying only these is sufficient.
#' @param by_time logical; should samples be ordered by acquisition time? 
#'   Ignored if \code{!is.null(md)} in which case samples will be ordered 
#'   as they are listed in \code{md[[md_cols$file]]}. (see details)
#'   
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#'   
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' 
#' @details 
#' By default, non-mass channels (e.g., time, event lengths) will be removed 
#' from the output SCE's assay data and instead stored in the object's internal 
#' cell metadata (\code{int_colData}) to assure these data are not subject to 
#' transformationsor other computations applied to the assay data.
#' 
#' For more than 1 sample, \code{prepData} will concatenate cells into a single 
#' \code{SingleCellExperiment} object. Note that cells will hereby be order by 
#' \code{"Time"}, regardless of whether \code{by_time = TRUE} or \code{FALSE}. 
#' Instead, \code{by_time} determines the sample (not cell!) order; 
#' i.e., whether samples should be kept in their original order, 
#' or should be re-ordered according to their acquision time 
#' stored in \code{keyword(flowSet, "$BTIM")}.
#' 
#' When a metadata table is specified (i.e. \code{!is.null(md)}), 
#' argument \code{by_time} will be ignored and sample ordering 
#' is instead determined by \code{md[[md_cols$file]]}.
#'   
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # channel-specific transformation
#' cf <- sample(seq_len(10)[-1], nrow(PBMC_panel), TRUE)
#' names(cf) <- PBMC_panel$fcs_colname
#' sce <- prepData(PBMC_fs, cofactor = cf)
#' int_metadata(sce)$cofactor
#' 
#' # input has different name for "condition"
#' md <- PBMC_md
#' m <- match("condition", names(md))
#' colnames(md)[m] <- "treatment"
#' 
#' # add additional factor variable batch ID
#' md$batch_id <- sample(c("A", "B"), nrow(md), TRUE)
#' 
#' # specify 'md_cols' that differ from defaults
#' factors <- list(factors = c("treatment", "batch_id"))
#' ei(prepData(PBMC_fs, PBMC_panel, md, md_cols = factors))
#' 
#' # without panel & metadata tables
#' sce <- prepData(raw_data)
#' 
#' # 'flowFrame' identifiers are used as sample IDs
#' levels(sce$sample_id)
#' 
#' # panel was guess with 'guessPanel';
#' # non-mass channels are set to marker class "none"
#' rowData(sce)
#' 
#' @importFrom methods new
#' @importFrom dplyr mutate_all rename
#' @importFrom flowCore colnames exprs exprs<- flowSet 
#'   fsApply identifier keyword parameters
#' @importFrom SingleCellExperiment SingleCellExperiment int_colData<- int_metadata<-
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom S4Vectors DataFrame
#' @export

prepData <- function(x, panel = NULL, md = NULL, 
    features = NULL, transform = TRUE, cofactor = 5,
    panel_cols = list(channel = "fcs_colname", 
        antigen = "antigen", class = "marker_class"),
    md_cols = list(
        file = "file_name", id = "sample_id", 
        factors = c("condition", "patient_id")),
    by_time = TRUE) {
    
    # read in data as 'flowSet'
    fs <- .read_fs(x)
    
    # reorder by acquisition time if 'md' is unspecified
    if (is.null(md) && by_time && length(fs) > 1) {
        ts <- keyword(fs, "$BTIM")
        if(any(vapply(ts, is.null, logical(1)))) {
            message("Not all samples contain information on their",
                " acquisition time; ignoring argument 'by_time'.",
                " Samples will be kept in their original order.")
        } else {
            fs <- fs[order(ts)]
        }
    }

    # assure panel & metadata are data.frames
    for (u in c("panel", "md"))
        if (!is.null(get(u)))
            assign(u, data.frame(get(u), 
                check.names = FALSE, 
                stringsAsFactors = FALSE))
    
    # fill up missing values with function defaults
    stopifnot(is.list(panel_cols), is.list(md_cols))
    args <- as.list(environment())
    for (i in c("md_cols", "panel_cols")) {
        defs <- as.list(formals("prepData")[[i]][-1])
        miss <- !names(defs) %in% names(args[[i]])    
        if (any(miss)) assign(i, c(args[[i]], defs[miss])[names(defs)])
    }
    stopifnot(
        c("channel", "antigen") %in% names(panel_cols),
        c("file", "id", "factors") %in% names(md_cols))
    
    # if unspecified, construct panel & metadata tables
    if (is.null(panel)) {
        panel <- guessPanel(fs[[1]])
        panel$marker_class <- ifelse(panel$use_channel, "state", "none")
    }
    if (is.null(md)) {
        ids <- fsApply(fs, identifier)
        md <- data.frame(
            file_name = ids,
            sample_id = basename(ids))
        md_cols$factors <- NULL
    }
    
    # check channels listed in panel 
    stopifnot(panel[[panel_cols$channel]] %in% colnames(fs))
    
    if (is.null(features)) {
        # default to channels listed in panel
        features <- as.character(panel[[panel_cols$channel]])
    } else {
        # check validity of 'features'
        chs <- colnames(fs)
        check1 <- is.logical(features) && length(features) == length(chs)
        check2 <- is.integer(features) && all(features %in% seq_along(chs))
        check3 <- all(features %in% chs)
        if (!any(check1, check2, check3))
            stop("Invalid argument 'features'. Should be either", 
                " a logial vector,\n  a numeric vector of indices, or",
                " a character vector of column names.")
        # subset panel & reorder features accordingly
        m <- match(panel[[panel_cols$channel]], features, nomatch = 0)
        panel <- panel[m != 0, , drop = FALSE]
        features <- features[m]
    }
    
    # check that filenames or identifiers 
    # match b/w 'flowSet' & metadata
    ids0 <- md[[md_cols$file]]
    ids1 <- fsApply(fs, identifier)
    ids2 <- keyword(fs, "FILENAME")
    if (length(unlist(ids2)) == length(fs))
        ids2 <- basename(ids2)
    check1 <- all(ids1 %in% ids0)
    check2 <- all(ids2 %in% ids0)
    ids_use <- which(c(check1, check2))[1]
    ids <- list(ids1, ids2)[[ids_use]]
    if (is.null(ids)) {
        stop("Couldn't match 'flowSet'/FCS filenames\n", 
            "with those listed in 'md[[md_cols$file]]'.")
    } else {
        # reorder 'flowSet' frames according to metadata table
        fs <- fs[match(md[[md_cols$file]], ids)]
    }
    
    # assure correctness of formats
    k <- c(md_cols$id, md_cols$factors)
    md <- md[, k, drop = FALSE] %>% 
        mutate_all(factor) %>% 
        rename("sample_id" = md_cols$id)
    
    # replace problematic characters
    as <- panel[[panel_cols$antigen]]
    as <- gsub("-", "_", as)
    as <- gsub(":", ".", as)
    as[is.na(as)] <- panel$fcs_colname[is.na(as)]
    
    # column & panel subsetting
    fs <- fs[, features]
    chs0 <- colnames(fs)
    
    # replace channel w/ antigen names
    m1 <- match(panel[[panel_cols$channel]], chs0, nomatch=0)
    m2 <- match(chs0, panel[[panel_cols$channel]], nomatch=0)
    as <- as[m2]
    ns <- table(as)
    for (a in names(ns)) if (ns[a] > 1)
        as[as == a] <- paste(a, seq_len(ns[a]), sep = ".")
    flowCore::colnames(fs)[m1] <- as
    chs <- colnames(fs)
    
    # get exprs.
    es <- matrix(fsApply(fs, exprs), byrow = TRUE,
        nrow = length(chs), dimnames = list(chs, NULL))
    
    # fix event times
    t <- grep("time", colnames(fs), ignore.case = TRUE)
    if (length(t) != 0) {
        ns <- fsApply(fs, nrow)
        t0 <- c(1, cumsum(ns) + 1)
        tx <- t0[-1] - 1
        for (i in seq_along(fs)[-1]) {
            idx <- seq(t0[i], tx[i])
            es[t, idx] <- es[t, idx] + es[t, tx[i - 1]]
        }
    }

    # get & check marker classes if provided
    mcs <- c("type", "state", "none")
    if (is.null(panel_cols$class) | is.null(panel[[panel_cols$class]])) {
        mcs <- factor("none", levels = mcs)
    } else {
        mcs <- factor(panel[[panel_cols$class]], levels = mcs)
        if (any(is.na(mcs)))
            stop("Invalid marker classes detected;",
                " valid classes are 'type', 'state', and 'none'.")
    }
    
    # construct row/colData & int_metadata
    rd <- DataFrame(
        row.names = chs, channel_name = chs0, 
        marker_name = chs, marker_class = mcs)
    m <- match(chs0, panel[[panel_cols$channel]], nomatch = 0)
    rd$use_channel <- panel$use_channel

    md$n_cells <- as.numeric(fsApply(fs, nrow))
    k <- setdiff(names(md), "n_cells")
    cd <- DataFrame(lapply(md[k], function(u) {
        v <- as.character(rep(u, md$n_cells))
        factor(v, levels = levels(u))
    }), row.names = NULL) 
    
    # construct SCE
    sce <- SingleCellExperiment(
        assays = list(counts = es), 
        rowData = rd, colData = cd,
        metadata = list(experiment_info = md))
    
    ds <- keyword(fs[[1]])
    keep <- lapply(c("\\$CYT$", "\\$CYTSN$"), grep, names(ds))
    int_metadata(sce)$description <- ds[unlist(keep)]
    
    # grep non-mass channels
    is_mass <- !is.na(.get_ms_from_chs(chs0))
    foo <- DataFrame(matrix(vector(), nrow = ncol(sce)))
    icd <- DataFrame(t(es[!is_mass, , drop = FALSE]), check.names = FALSE)
    colnames(icd) <- rownames(es)[!is_mass]
    icd$reducedDims <- icd$altExps <- foo
    # store & exclude from assay data
    int_colData(sce) <- icd
    sce <- sce[is_mass, ]
    
    # (optionally) do arcsinh-transformation & return SCE
    if (transform) .transform(sce, cofactor) else sce
}
