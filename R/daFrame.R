# ==============================================================================
# daFrame constructor
# ------------------------------------------------------------------------------
#' @rdname daFrame-class
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
#' @param cols_to_use a logical vector, numeric vector of column indices,
#'   or character vector of channel names. Specified which column to keep 
#'   from the input data. Defaults to the channels listed in the input panel.
#' @param cofactor numeric cofactor to use for arcsinh-transformation.
#' @param panel_cols a named list specifying the column names of \code{md}
#'   that contain the FCS file names, sample IDs, and factors of interest
#'   (batch, condition, treatment etc.).
#' @param md_cols a names list specifying the column names of \code{panel}
#'   that contain the channel names, targeted protein markers, and (optionally) 
#'   marker classes. 
#'   
#' @return an object of class \code{daFrame}.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' @importFrom dplyr %>% mutate_at
#' @importFrom flowCore colnames exprs exprs<- flowSet 
#'   fsApply identifier isFCSfile keyword read.flowSet
#' @importFrom SingleCellExperiment SingleCellExperiment

setMethod("daFrame",
    signature(x="flowSet", panel="data.frame", md="data.frame"),
    function(x, panel, md, cols_to_use=NULL, cofactor=5,
        panel_cols=list(
            channel="fcs_colname", antigen="antigen", class="marker_class"),
        md_cols=list(
            file="file_name", id="sample_id", 
            factors=c("condition", "patient_id"))) {

        # check validity of input arguments
        stopifnot(is.list(panel_cols), is.list(md_cols),
            all(c("channel", "antigen", "class") %in% names(panel_cols)),
            all(c("file", "id", "factors") %in% names(md_cols)))
        stopifnot(is.numeric(cofactor), length(cofactor) == 1, cofactor > 0)
        fs <- x
        
        # check channels listed in panel 
        stopifnot(all(panel[[panel_cols$channel]] %in% colnames(fs)))
        
        if (is.null(cols_to_use)) {
            # default to channels listed in panel
            cols_to_use <- as.character(panel[[panel_cols$channel]])
        } else {
            # check validity of 'cols_to_use'
            chs <- colnames(fs)
            check1 <- is.logical(cols_to_use) && 
                length(cols_to_use) == length(chs)
            check2 <- all(cols_to_use %in% chs)
            check3 <- is.integer(cols_to_use) && 
                all(cols_to_use %in% seq_along(chs))
            if (!any(check1, check2, check3))
                stop("Invalid argument 'cols_to_use'. Should be either", 
                    " a logial vector,\n  a numeric vector of indices, or",
                    " a character vector of column names.")
        }
        
        # use identifiers if filenames are not specified
        ids <- c(keyword(fs, "FILENAME"))
        if (is.null(unlist(ids)))
            ids <- c(fsApply(fs, identifier))
        
        # check that filenames match b/w flowSet & metadata &
        # reorder flowSet according to metadata table
        stopifnot(all(ids %in% md[[md_cols$file]]))
        fs <- fs[match(ids, md[[md_cols$file]])]
        
        # transformation
        fs <- fsApply(fs, function(ff) {
            exprs(ff) <- asinh(exprs(ff) / cofactor)
            return(ff)
        })
        
        # assure correctness of formats
        md <- data.frame(md)
        md <- md %>% mutate_at(c(md_cols$id, md_cols$factors), factor)
        
        # replace problematic characters
        antigens <- panel[[panel_cols$antigen]]
        antigens <- gsub("-", "_", antigens)
        antigens <- gsub(":", ".", antigens)
        
        # column subsetting
        fs <- fs[, cols_to_use]
        chs0 <- colnames(fs)
        
        # replace channel w/ antigen names
        m1 <- match(panel[[panel_cols$channel]], chs0, nomatch=0)
        m2 <- match(chs0, panel[[panel_cols$channel]], nomatch=0)
        flowCore::colnames(fs)[m1] <- antigens[m2]
        chs <- colnames(fs)
        
        # get exprs.
        es <- matrix(fsApply(fs, exprs), 
            ncol=length(chs), dimnames=list(NULL, chs))
        
        # get nb. of cells per sample
        n_cells <- as.numeric(fsApply(fs, nrow))
        n_cells <- setNames(n_cells, md[[md_cols$id]])
        
        # get & check marker classes if provided
        valid_mcs <- c("type", "state", "none")
        if (is.null(panel[[panel_cols$class]])) {
            mcs <- factor("none", levels=valid_mcs)
        } else {
            mcs <- factor(panel[[panel_cols$class]], levels=valid_mcs)
            mcs <- mcs[match(chs0, panel[[panel_cols$channel]])]
            if (any(is.na(mcs)))
                stop("Invalid marker classes detected.",
                    " Valid classes are 'type', 'state', and 'none'.")
        }
        
        # construct row & column data
        k <- setdiff(names(md), md_cols$file)
        row_data <- lapply(md[k], function(u) {
            v <- as.character(rep(u, n_cells))
            factor(v, levels = levels(u))
        }) %>% data.frame(row.names = NULL)
        
        col_data <- data.frame(
            row.names=chs, channel_name=chs0, 
            marker_name=chs, marker_class=mcs)
        
        # construct daFrame
        sce <- SingleCellExperiment(
            assays=list(exprs=es), 
            rowData=row_data, 
            colData=col_data,
            metadata=list(
                experiment_info=md, 
                n_cells=n_cells, 
                cofactor=cofactor))
        new("daFrame", sce)
    }
)
# ------------------------------------------------------------------------------

#' @rdname daFrame-class
setMethod("daFrame",
    signature(x="character", panel="data.frame", md="data.frame"),
    function(x, panel, md, ...) {
        stopifnot(dir.exists(x))
        fcs <- list.files(x, ".fcs$", full.names=TRUE, ignore.case=TRUE)
        if (length(fcs) == 1) 
            stop("The specified directory contains only a single FCS file.")
        stopifnot(all(vapply(fcs, isFCSfile, logical(1))))
        fs <- read.flowSet(fcs, transformation=FALSE, truncate_max_range=FALSE)
        for (i in seq_along(fs))
            description(fs[[i]])$FILENAME <- description(fs[[i]])$GUID
        daFrame(fs, panel, md, ...)
    }
)
# ------------------------------------------------------------------------------

#' @rdname daFrame-class
setMethod("daFrame",
    signature(x="daFrame", panel="matrix", md="ANY"),
    function(x, panel, md, ...) daFrame(x, data.frame(panel), md, ...))

#' @rdname daFrame-class
setMethod("daFrame",
    signature(x="daFrame", panel="ANY", md="matrix"),
    function(x, panel, md, ...) daFrame(x, panel, data.frame(md), ...))