# ==============================================================================
# Debarcoding frame class
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
#' \code{\link{assignPrelim}} will return a \code{dbFrame} containing the
#' input measurement data, barcoding scheme, and preliminary event assignments.
#' \item assignments will be made final by \code{\link{applyCutoffs}}.
#' Optionally, population-specific separation cutoffs may be estimated 
#' by running \code{\link{estCutoffs}} prior to this.
#' \item \code{\link{plotYields}}, \code{\link{plotEvents}} and 
#' \code{\link{plotMahal}} aim to guide selection of devoncolution parameters 
#' and to give a sense of the resulting barcode assignment quality.
#' }
#' \code{show(dbFrame)} will display \itemize{
#' \item the dimensionality of the measurement data and number of barcodes
#' \item current assignments in order of decreasing population size
#' \item current separation cutoffs
#' \item the average and per-population yield 
#'       that will be achieven upon debarcoding}
#' 
#' @slot exprs  
#' a matrix containing raw intensities of the input flowFrame.
#' @slot bc_key 
#' binary barcoding scheme with numeric masses as column names 
#' and samples names as row names OR a numeric vector of barcode masses.
#' @slot bc_ids
#' vector of barcode IDs. If a barcoding scheme is supplied, 
#' the respective binary code's row name, else, the mass of the respective 
#' barcode channel.
#' @slot deltas 
#' numeric vector of separations between positive and negative 
#' barcode populations computed from normalized barcode intensities.
#' @slot normed_bcs 
#' matrix containing normalized barcode intensities.
#' @slot mhl_dists
#' mahalanobis distances.
#' @slot sep_cutoffs
#' numeric vector of distance separation cutoffs between positive and negative 
#' barcode populations above which events will be unassigned.
#' @slot mhl_cutoff
#' non-negative and non-zero numeric value specifying the Mahalanobis distance 
#' below which events will be unassigned.
#' @slot counts
#' matrix of dimension (# barcodes)x(101) where each row contains the number 
#' of events within a barcode for which positive and negative populations 
#' are separated by a distance between in [0,0.01), ..., [0.99,1], respectively.
#' @slot yields
#' a matrix of dimension (# barcodes)x(101) where each row contains the 
#' percentage of events within a barcode that will be obtained after applying
#' a separation cutoff of 0, 0.01, ..., 1, respectively.
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom methods new
#' @export
# ------------------------------------------------------------------------------

dbFrame <- setClass(Class="dbFrame", package="CATALYST", slots=c(
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

# ------------------------------------------------------------------------------

setValidity(Class="dbFrame", 
    method=function(object){
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
    })

# ==============================================================================
# Differential analysis frame class
# ------------------------------------------------------------------------------
#' @rdname daFrame-class
#' @name daFrame-class
#' 
#' @title Differential analysis frame class
#' @description 
#' Represents the data returned by and used throughout differential analysis.
#' 
#' @slot assays 
#' list of length one containing the arcsinh-transformed expressions.
#' @slot rowData 
#' the metadata information for each event, and its cluster ID
#' as inferred by the initial \code{\link{FlowSOM}} clustering.
#' @slot colData 
#' a data.frame with the following columns:\itemize{
#' \item \code{channel} original column name in the input \code{flowSet}
#' \item \code{type1}, \code{type2} logical vectors indicating, 
#' for each antigen, whether it was used for clustering}
#' @slot metadata 
#' a named list containing:\itemize{
#' \item the original metadata-table
#' \item \code{panel}: the original panel-table
#' \item \code{n_events}: the number of events measured per sample
#' \item \code{SOM_codes}: a 100 x p matrix of SOM codes, 
#' where p = no. of measurement parameters
#' \item \code{cluster_codes}: cluster codes for the initial 
#' \code{\link{FlowSOM}} clustering (column \code{'100'}), the 
#' \code{\link{ConsensusClusterPlus}} metaclustering (columns \code{'2'} 
#' through {'20'}, and mergings as labeled by \code{\link{mergeClusters}}}
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ConsensusClusterPlus Rtsne SummarizedExperiment
#' @importFrom flowCore colnames exprs fsApply parameters pData
#' @importFrom FlowSOM BuildSOM ReadInput
#' @importFrom methods new
#' @importFrom S4Vectors DataFrame SimpleList
#' @export

# ------------------------------------------------------------------------------
# class definition
setClass(
    Class="daFrame", 
    contains="SummarizedExperiment")

# ------------------------------------------------------------------------------
# constructor
#' @rdname daFrame-class
#' 
#' @param fs a \code{\link{flowSet}} holding all samples.
#' @param panel a data.frame with the following columns:
#' \itemize{
#' \item \code{fcs_colname}: the marker's column name in the FCS file
#' \item \code{antigen}: the targeted protein marker}
#' @param md a data.frame with columns describing the experiment (experimental condition, batch, sample identifiers, etc.)
#' \itemize{
#' \item \code{file_name}: the FCS file name
#' \item \code{sample_id}: a unique sample identifier
#' \item \code{condition}: brief sample description (e.g. REF)
#' \item \code{patient_id}: the patient ID}
#' 
#' 
#' @param cols_to_use a logical vector OR numeric vector of indices OR character vector of column names.
#' Specifies the columns to keep from the input \code{flowSet}.
#' @param cofactor cofactor to use for arcsinh-transformation.
#' 
#' @export
#' @import SummarizedExperiment

daFrame <- function(fs, panel, md, cols_to_use=NULL, cofactor=5, 
                    exprs_colname="fcs_colname",
                    antigen_colname="antigen") {
    
    # set/check colnames of panel 
    chs <- flowCore::colnames(fs)
    if (is.null(cols_to_use))
        cols_to_use <- chs
    check_validity_columns(cols_to_use, chs)
    
    # replace problematic characters
    check_validity_columns(exprs_colname, colnames(panel))
    check_validity_columns(antigen_colname, colnames(panel))
    antigens <- gsub("-", "_", panel[[antigen_colname]])
    
    # arcsinh-transformation & column subsetting
    fs <- fs[, cols_to_use]
    fs <- fsApply(fs, function(ff) {
        flowCore::exprs(ff) <- asinh(exprs(ff)/cofactor)
        ff
    })
    
    md <- data.frame(md)
    chs <- flowCore::colnames(fs)
    m1 <- match(panel[[exprs_colname]], chs, nomatch=0)
    m2 <- match(chs, panel[[exprs_colname]])
    flowCore::colnames(fs)[m1] <- antigens[m2]
    es <- matrix(fsApply(fs, exprs), 
                 ncol=length(chs),
                 dimnames=list(NULL, flowCore::colnames(fs)))
    n_events <- fsApply(fs, nrow)
    n_events <- setNames(as.numeric(n_events), md$sample_id)

    # construct SummarizedExperiment
    conditions <- grep("condition", colnames(md), value=TRUE)
    conditions <- sapply(conditions, function(i) rep(md[, i], n_events))
    row_data <- S4Vectors::DataFrame(sample_id=rep(md$sample_id, n_events), conditions)
    col_data <- S4Vectors::DataFrame(channel=chs, row.names=colnames(es))

    new("daFrame", 
        SummarizedExperiment(assays=SimpleList(es=es),
                             rowData=row_data, 
                             colData=col_data,
                             metadata=list(design=md, n_events=n_events)))
}

# validity
# ------------------------------------------------------------------------------
# setValidity(Class="daFrame", 
#     method=function(object){
#         # check panel & metadata column names
#         if (5 != sum(colnames(panel(object)) %in% 
#                 c("Metal", "Isotope", "Antigen", "Lineage", "Functional")))
#             return(message(""))
#         if (3 != sum(colnames(metadata(object)) %in% 
#                 c("file_name", "sample_id", "condition")))
#             return(message(""))
#         # check that all flowSet file names 
#         # occur in metadata file_name column
#         if (!all(keyword(fs, "FILENAME") %in% md$file_name))
#             return(message(""))
#         return(TRUE)
#     })
