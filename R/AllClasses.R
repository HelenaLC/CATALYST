# ==============================================================================
# Debarcoding frame class definition
# ------------------------------------------------------------------------------
#' @rdname dbFrame-class
#' @name dbFrame-class
#' 
#' @title Debarcoding frame class
#' @description 
#' This class represents the data returned by and used throughout debarcoding.
#' 
#' @details 
#' Objects of class \code{dbFrame} hold all data required for debarcoding:
#' \enumerate{
#' \item as the initial step of single-cell deconcolution, 
#'   \code{\link{assignPrelim}} will return a \code{dbFrame} containing
#'   the input measurement data, barcoding scheme, and preliminary assignments.
#' \item assignments will be made final by \code{\link{applyCutoffs}}.
#'   Optionally, population-specific separation cutoffs may be estimated 
#'   by running \code{\link{estCutoffs}} prior to this.
#' \item \code{\link{plotYields}}, \code{\link{plotEvents}} and 
#'   \code{\link{plotMahal}} aim to guide devoncolution parameter selection,
#'   and to give a sense of the resulting barcode assignment quality.}
#' \code{show(dbFrame)} will display \itemize{
#' \item the dimensionality of the measurement data and number of barcodes
#' \item current assignments in order of decreasing population size
#' \item current separation cutoffs
#' \item the mean & per-population yield that'll be achieved upon debarcoding}
#' 
#' @slot exprs  
#'   a matrix containing raw intensities of the input flowFrame.
#' @slot bc_key 
#'   binary barcoding scheme with numeric masses as column names and
#'   samples names as row names OR a numeric vector of barcode masses.
#' @slot bc_ids
#'   vector of barcode IDs. If a barcoding scheme is supplied, the respective 
#'   binary code's row name, else, the mass of the respective barcode channel.
#' @slot deltas 
#'   numeric vector of separations between positive and negative 
#'   barcode populations computed from normalized barcode intensities.
#' @slot normed_bcs 
#'   matrix containing normalized barcode intensities.
#' @slot mhl_dists
#'   mahalanobis distances.
#' @slot sep_cutoffs
#'   numeric vector of distance separation cutoffs between positive 
#'   and negative barcode populations above which events will be unassigned.
#' @slot mhl_cutoff
#'   non-negative and non-zero numeric value specifying the 
#'   Mahalanobis distance below which events will be unassigned.
#' @slot counts
#'   matrix of dimension (# barcodes)x(101) where each row contains the number 
#'   of events within a barcode for which positive & negative populations are 
#'   separated by a distance between in [0,0.01), ..., [0.99,1], respectively.
#' @slot yields
#'   a matrix of dimension (# barcodes)x(101) where each row contains the 
#'   percentage of events within a barcode that will be obtained after applying
#'   a separation cutoff of 0, 0.01, ..., 1, respectively.
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom methods new
#' @export

dbFrame <- setClass(
    Class="dbFrame", 
    package="CATALYST", 
    slots=c(
        exprs="matrix",
        bc_key="data.frame",
        bc_ids="vector",
        deltas="numeric",
        normed_bcs ="matrix",
        mhl_dists = "numeric",
        sep_cutoffs="numeric",
        mhl_cutoff="numeric",
        counts="matrix",
        yields="matrix"))

# validity check
setValidity(Class="dbFrame", 
    method=function(object) {
        n <- nrow(exprs(object))
        ms <- gsub("[[:alpha:][:punct:]]", "", colnames(exprs(object)))
        # check that all barcode masses occur in the measurement data
        if (!all(colnames(bc_key(object)) %in% ms))
            return(message("Invalid 'bc_key': Column names must be numeric",
                "\nand coherent with masses eobjecttracted from 'exprs'."))
        # check that 'normed_bcs' is of dimension
        # number events x number barcode channels
        if (!all.equal(dim(normed_bcs(object)), c(n, ncol(bc_key(object)))))
            return(message("'normed_bcs' should be of dimension",
                "nrow('exprs') x ncol('bc_key')."))
        # check that 'bc_ids', 'deltas' and 'mhl_dists' 
        # are of length number of events
        if (!all(c(length(bc_ids(object)), length(deltas(object))) == n))
            return(message("'bc_ids' and 'deltas' should have\n",
                "as many entries as numbers of rows in 'exprs'."))
        if (!length(mhl_dists(object)) %in% c(0, n))
            return(message("'mhl_dists' should have\n",
                "as many entries as numbers of rows in 'exprs'."))
        # check that all 'bc_ids" are 0 = "unassigned"
        # or occur as row names in the 'bc_key'
        if ((valid <- sum(bc_ids(object) %in% 
                c(0, rownames(bc_key(object))))) != n)
            return(message(n-valid, "/", n, " 'bc_ids' are invalid.\n",
                "'bc_ids' should be either 0 = \"unassigned\"\n",
                "or occur as rownames in the 'bc_key'."))
        return(TRUE)
    }
)

# ==============================================================================
# Differential analysis frame class definition
# ------------------------------------------------------------------------------
#' @rdname daFrame-class
#' @name daFrame-class
#' @title Differential analysis frame class
#' 
#' @description Represents the data returned by and used throughout differential analysis.
#' 
#' \describe{
#' \item{\code{assays}}{a list of length 1 containing the measurement data.}
#' \item{\code{rowData}}{a \code{\link{DataFrame}} containing the row metadata.}
#' \item{\code{colData}}{a \code{\link{DataFrame}} containing the column metadata.}
#' \item{\code{metadata}}{a named list containing: \itemize{
#'   \item{\code{experimental_design}: the original metadata table.}
#'   \item{\code{n_cells}: number of events (cells) measured per sample.}
#'   \item{\code{cofactor}: numeric cofactor used for arcsinh-transformation.}
#'   \item{\code{cluster_codes}:
#'     cluster codes for the initial \pkg{FlowSOM} clustering (\code{"som100"}), 
#'     the \pkg{ConsensusClusterPlus} metaclustering (\code{"metaX"}), and all
#'     manual mergings done with \code{\link{mergeClusters}}.}
#'   \item{\code{delta_area}: the Delta Area plot 
#'     (see \pkg{ConsensusClusterPlus} for details).}
#'   \item{\code{experimental_design}: the original metadata table}}}}
#'       
#' @author Helena Lucia Crowell \email{helena@crowells.eu}
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods new
#' 
#' @export

setClass(Class="daFrame", package="CATALYST", contains="SummarizedExperiment")

# ==============================================================================
# daFrame constructor
# ------------------------------------------------------------------------------
#' @rdname daFrame-class
#' 
#' @param x a \code{flowSet} holding all samples or a path to a set of FCS files.
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
#' @importFrom flowCore colnames exprs exprs<- fsApply isFCSfile keyword read.flowSet
#' @importFrom SummarizedExperiment SummarizedExperiment
#' 
#' @export

daFrame <- function(x, panel, md, cols_to_use=NULL, cofactor=5,
    panel_cols=list(channel="fcs_colname", antigen="antigen", class="marker_class"),
    md_cols=list(file="file_name", id="sample_id", factors=c("condition", "patient_id"))) {
    
    # check validity of input arguments
    stopifnot(is.numeric(cofactor), length(cofactor) == 1, cofactor > 0)
    stopifnot(is.list(panel_cols), is.list(md_cols),
        all(c("channel", "antigen", "class") %in% names(panel_cols)),
        all(c("file", "id", "factors") %in% names(md_cols)))
    
    # get input data
    if (is.character(x)) {
        stopifnot(dir.exists(x))
        fcs <- list.files(x, ".fcs$", full.names=TRUE, ignore.case=TRUE)
        if (length(fcs) == 1) 
            stop("The specified directory contains only a single FCS file.")
        stopifnot(all(vapply(fcs, flowCore::isFCSfile, logical(1))))
        fs <- flowCore::read.flowSet(fcs, transformation=FALSE, truncate_max_range=FALSE)
    } else if (class(x) == "flowSet") {
        fs <- x
    } else {
        stop("'x' should be either a flowSet or\na character string",
            " that specifyies a path to a set of FCS files.")
    }
    
    # check channels listed in panel 
    chs <- flowCore::colnames(fs)
    stopifnot(all(panel[[panel_cols$channel]] %in% chs))
    
    # default to channels listed in panel
    if (is.null(cols_to_use)) {
        cols_to_use <- as.character(panel[[panel_cols$channel]])
    } else {
        # check validity of 'cols_to_use'
        check1 <- is.logical(cols_to_use) && length(cols_to_use) == length(chs)
        check2 <- all(cols_to_use %in% chs)
        check3 <- is.integer(cols_to_use) && all(cols_to_use %in% seq_along(chs))
        if (!any(check1, check2, check3))
            stop("Invalid argument 'cols_to_use'. Should be either", 
                " a logial vector,\n  a numeric vector of indices, or",
                " a character vector of column names.")
    }
    
    # check that filenames match b/w flowSet & metadata
    stopifnot(all(flowCore::keyword(fs, "FILENAME") %in% md[[md_cols$file]]))
    
    # reorder flowSet according to metadata table & 
    m <- match(flowCore::keyword(fs, "FILENAME"), md[[md_cols$file]])
    fs <- fs[m]
    
    # transformation
    fs <- flowCore::fsApply(fs, function(ff) {
        flowCore::exprs(ff) <- asinh(flowCore::exprs(ff) / cofactor)
        return(ff)
    })
    
    # assure correctness of formats
    md <- data.frame(md)
    for (i in md_cols$factors) 
        md[[i]] <- factor(md[[i]])
    
    # replace problematic characters
    antigens <- panel[[panel_cols$antigen]]
    antigens <- gsub("-", "_", antigens)
    antigens <- gsub(":", ".", antigens)
    
    # column subsetting
    fs <- fs[, cols_to_use]
    chs <- flowCore::colnames(fs)
    
    # replace channel w/ antigen names
    m1 <- match(panel[[panel_cols$channel]], chs, nomatch=0)
    m2 <- match(chs, panel[[panel_cols$channel]], nomatch=0)
    flowCore::colnames(fs)[m1] <- antigens[m2]
    chs <- flowCore::colnames(fs)
    
    # get exprs.
    es <- matrix(flowCore::fsApply(fs, exprs), 
        ncol=length(chs), dimnames=list(NULL, chs))
    
    # get nb. of cells per sample
    n_cells <- flowCore::fsApply(fs, nrow)
    n_cells <- setNames(as.numeric(n_cells), md[[md_cols$id]])
    
    # get & check marker classes if provided
    valid_mcs <- c("type", "state", "none")
    if (is.null(panel[[panel_cols$class]])) {
        mcs <- factor("none", levels=valid_mcs)
    } else {
        mcs <- factor(panel[[panel_cols$class]])
        if (!all(mcs %in% valid_mcs))
            stop("Invalid marker classes detected.",
                " Valid classes are 'type', 'state', and 'none'.")
        levels(mcs) <- valid_mcs
    }
    
    # construct row & column data
    row_data <- data.frame(row.names=NULL,
        apply(md, 2, rep, n_cells)[, names(md) != md_cols$file])
    col_data <- data.frame(row.names=chs, 
        channel_name=chs, marker_name=chs, marker_class=mcs)
    
    # construct daFrame
    new("daFrame", SummarizedExperiment(
        assays=list(exprs=es), rowData=row_data, colData=col_data,
        metadata=list(experiment_info=md, n_cells=n_cells, cofactor=cofactor)))
}

# validity check
setValidity(Class="daFrame", 
    method=function(object) {
        x <- deparse(substitute(object))
        # ----------------------------------------------------------------------
        # check colData(x)$marker_class
        valid_classes <- c("type", "state", "none")
        lvls <- levels(colData(object)$marker_class)
        if (any(is.na(match(lvls, valid_classes)))
            | any(!colData(object)$marker_class %in% lvls))
            return(message("colData(", x, ")$marker_class ",
                "should be of type factor\nwith levels ", 
                paste(dQuote(valid_classes), collapse=", "), "."))
        # ----------------------------------------------------------------------
        # check names(metadata(x))
        md_nms <- c("experiment_info", "n_cells", "cofactor",
            "SOM_codes", "cluster_codes", "delta_area")
        if (!all(names(metadata(object)) %in% md_nms))
            return(message("Invalid names(metadata(", x, ").", 
                " Metadata should contain:\n", 
                paste(dQuote(md_nms), collapse=", "), "."))
        # ----------------------------------------------------------------------
        # check metadata(x)$n_cells
        if (nrow(object) != sum(metadata(object)$n_cells))
            return(message("nrow(", x, ") != sum(metadata(", x, ")$n_cells)"))
    }
)
