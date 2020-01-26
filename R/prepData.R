#' @rdname prepData
#' @title Data preparation
#' 
#' @param x a \code{flowSet} holding all samples 
#'   or a path to a set of FCS files.
#' @param panel a data.frame containing, for each channel, 
#'   its column name in the input data, targeted protein marker,
#'   and (optionally) class ("type", "state", or "none").
#' @param md a table with column describing the experiment.
#'   An exemplary metadata table could look as follows: \itemize{
#'   \item\code{file_name}: the FCS file name
#'   \item\code{sample_id}: a unique sample identifier
#'   \item\code{patient_id}: the patient ID
#'   \item\code{condition}: brief sample description 
#'     (e.g. reference/stimulated, healthy/diseased)}
#' @param features a logical vector, numeric vector of column indices,
#'   or character vector of channel names. Specified which column to keep 
#'   from the input data. Defaults to the channels listed in the input panel.
#' @param cofactor numeric cofactor to use for arcsinh-transformation.
#' @param panel_cols a names list specifying the column names of \code{panel}
#'   that contain the channel names, targeted protein markers, and (optionally) 
#'   marker classes. 
#' @param md_cols a named list specifying the column names of \code{md}
#'   that contain the FCS file names, sample IDs, and factors of interest
#'   (batch, condition, treatment etc.).
#'   
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#'   
#' @return an object of class \code{daFrame}.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' prepData(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' @importFrom methods new
#' @importFrom dplyr mutate_all rename
#' @importFrom flowCore colnames exprs exprs<- flowSet 
#'   fsApply identifier isFCSfile keyword read.flowSet
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @export

prepData <- function(x, panel, md, features=NULL, cofactor=5,
    panel_cols=list(channel="fcs_colname", antigen="antigen", class="marker_class"),
    md_cols=list(file="file_name", id="sample_id", factors=c("condition", "patient_id"))) {
    
    # check validity of input arguments
    for (u in c("panel", "md"))
        if (!is(v <- get(u), "data.frame"))
            assign(u, data.frame(v, 
                check.names = FALSE, 
                stringsAsFactors = FALSE))
    
    stopifnot(is.list(panel_cols), is.list(md_cols),
        c("channel", "antigen") %in% names(panel_cols),
        c("file", "id", "factors") %in% names(md_cols))

    if (is(x, "flowSet")) {
        fs <- x
    } else if (is.character(x)) {
        stopifnot(dir.exists(x))
        fcs <- list.files(x, ".fcs$", full.names=TRUE, ignore.case=TRUE)
        if (length(fcs) < 2) 
            stop("The specified directory contains",
                " none or only a single FCS file.")
        stopifnot(all(vapply(fcs, isFCSfile, logical(1))))
        fs <- read.flowSet(fcs, transformation=FALSE, truncate_max_range=FALSE)
    } else {
        stop("Invalid argument 'x'; should be either a flowSet",
            " or a character string specifying the path to", 
            " a directory containing a set of FCS files.")
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
        m <- match(panel[[panel_cols$channel]], features, nomatch = 0)
        features <- features[m]
    }
    
    if (!is.null(cofactor))
        stopifnot(is.numeric(cofactor), cofactor > 0,
            length(cofactor) %in% c(1, length(features)))
    
    # use identifiers if filenames are not specified
    ids <- c(keyword(fs, "FILENAME"))
    if (is.null(unlist(ids))) 
        ids <- c(fsApply(fs, identifier))
    
    # check that filenames match b/w flowSet & metadata &
    # reorder flowSet according to metadata table
    stopifnot(all(ids %in% md[[md_cols$file]]))
    fs <- fs[m <- match(ids, md[[md_cols$file]])]
    
    # transformation
    if (!is.null(cofactor)) {
        if (length(cofactor) == 1) {
            cofactor <- rep(cofactor, ncol(fs[[1]]))
        } else if (!is.null(names(cofactor))) {
            stopifnot(features %in% names(cofactor))
            m <- match(features, names(cofactor))
            cofactor <- cofactor[m]
        } 
        fs <- fsApply(fs, function(ff) {
            exprs(ff) <- asinh(sweep(exprs(ff), 2, cofactor, "/"))
            #exprs(ff) <- asinh(exprs(ff) / cofactor)
            return(ff)
        })
    }
    
    # assure correctness of formats
    k <- c(md_cols$id, md_cols$factors)
    md <- data.frame(md)[, k] %>% 
        mutate_all(factor) %>% 
        dplyr::rename("sample_id" = md_cols$id)
    o <- order(md[[md_cols$factors[1]]])
    md$sample_id <- factor(md$sample_id, levels = md$sample_id[o])
    
    # replace problematic characters
    antigens <- panel[[panel_cols$antigen]]
    antigens <- gsub("-", "_", antigens)
    antigens <- gsub(":", ".", antigens)
    
    # column subsetting
    fs <- fs[, features]
    chs0 <- colnames(fs)
    
    # replace channel w/ antigen names
    m1 <- match(panel[[panel_cols$channel]], chs0, nomatch=0)
    m2 <- match(chs0, panel[[panel_cols$channel]], nomatch=0)
    flowCore::colnames(fs)[m1] <- antigens[m2]
    chs <- colnames(fs)
    
    # get exprs.
    es <- matrix(fsApply(fs, exprs), byrow = TRUE,
        nrow=length(chs), dimnames=list(chs, NULL))
    # get nb. of cells per sample
    md$n_cells <- as.numeric(fsApply(fs, nrow))
    
    # get & check marker classes if provided
    valid_mcs <- c("type", "state", "none")
    if (is.null(panel_cols$class) | is.null(panel[[panel_cols$class]])) {
        mcs <- factor("none", levels=valid_mcs)
    } else {
        mcs <- factor(panel[[panel_cols$class]], levels=valid_mcs)
        mcs <- mcs[match(chs0, panel[[panel_cols$channel]])]
        if (any(is.na(mcs)))
            stop("Invalid marker classes detected.",
                " Valid classes are 'type', 'state', and 'none'.")
    }
    
    # construct row & column data
    rd <- DataFrame(
        row.names=chs, channel_name=chs0, 
        marker_name=chs, marker_class=mcs)
    
    k <- setdiff(names(md), "n_cells")
    cd <- DataFrame(lapply(md[k], function(u) {
        v <- as.character(rep(u, md$n_cells))
        factor(v, levels = levels(u))
    }), row.names = NULL) 
    
    # construct SCE
    SingleCellExperiment(
        assays=list(exprs=es), rowData=rd, colData=cd,
        metadata=list(experiment_info=md, cofactor=cofactor))
}
