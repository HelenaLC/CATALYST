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
#' list of length one containing cofactor 5 arcsinh-transformed expressions 
#' of lineage and functional markers as specified in the metadata-table.
#' @slot rowData 
#' data.frame containing cluster codes (column \code{cluster_codes}) and IDs 
#' (column \code{cluster_ids}) for the initial \code{\link{FlowSOM}}-clustering 
#' (k=100), and \code{\link{ConsensusClusterPlus}} metaclustering (k=20-2).
#' @slot colData 
#' binary table indicating, for each antigen, 
#' whether it is a lineage or functional marker.
#' @slot metadata 
#' contains the original metadata- and panel-table, the number of events 
#' measured per sample (\code{n_events}), a 100 x p matrix of SOM codes, 
#' where n = no. of lineage + no. of functional markers
#' 
#' @details Objects of class \code{daFrame} 
#' hold all data required for differential analysis:
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ConsensusClusterPlus Rtsne SummarizedExperiment
#' @importFrom flowCore colnames exprs fsApply parameters pData
#' @importFrom FlowSOM BuildSOM ReadInput
#' @importFrom methods new
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
#' \item \code{file_name}: the .fcs file name
#' \item \code{sample_id}: a unique sample identifier
#' \item \code{condition}: brief sample description (e.g. REF)
#' \item \code{patient_id}: the patient ID}
#' 
#' @param md a data.frame with the following columns:
#' \itemize{
#' \item \code{Metal} and \code{Isotope}: 
#' symbol and atomic mass of the element conjugated to the antibody
#' \item \code{Antigen}: the targeted protein marker
#' \item \code{Lineage} and \code{Functional}: 0 or 1 to indicate 
#' whether an antibody is a lineage or funcitonal marker}
#' 
#' @export

daFrame <- function(fs, panel, md) {
    # replace problematic characters
    pd <- pData(parameters(fs[[1]]))
    pd$desc <- gsub("-", "_", pd$desc)
    panel$Antigen <- gsub("-", "_", panel$Antigen)
    
    # arcsinh-transformation & column subsetting
    l <- panel$Antigen[as.logical(panel$Lineage)]
    f <- panel$Antigen[as.logical(panel$Functional)]
    fs <- fsApply(fs, function(ff) {
        flowCore::colnames(ff) <- pd$desc
        flowCore::exprs(ff) <- asinh(exprs(ff[, c(l,f)])/5)
        return(ff)
    })
    es <- fsApply(fs, exprs)
    n_events <- fsApply(fs, nrow)
    n_events <- setNames(as.numeric(n_events), md$sample_id)
    
    # flowSOM clustering
    message("o running FlowSOM clustering...")
    fsom <- ReadInput(fs, transform=FALSE, scale=FALSE)
    som <- BuildSOM(fsom, colsToUse=l, silent=TRUE)
    codes <- som$map$codes
    cluster_ids <- som$map$mapping[, 1]
    
    # metaclustering
    message("o running ConsensusClusterPlus metaclustering...")
    pdf(NULL)
    mc <- suppressMessages(ConsensusClusterPlus(t(codes), 
        maxK=20, reps=100, distance="euclidean", plot="pdf"))
    dev.off()
    
    # get cluster codes for k = 100, 2-20
    cluster_codes <- data.frame(matrix(0, 100, 20, 
        dimnames=list(NULL, c(100, 2:20))), check.names=FALSE)
    for (k in seq_len(20)[-1])
        cluster_codes[, k] <- mc[[k]]$consensusClass
    cluster_codes$`100` <- seq_len(100)
    
    # construct SummarizedExperiment
    inds <- panel$Lineage | panel$Functional
    row_data <- data.frame(
        condition=rep(md$condition, n_events),
        sample_id=rep(md$sample_id, n_events),
        cluster_id=cluster_ids)
    col_data <- data.frame(
        lineage=panel$Lineage[inds],
        functional=panel$Functional[inds],
        row.names=colnames(es))
    metadata <- list(
        md, 
        panel=panel, 
        n_events=n_events,
        SOM_codes=codes,
        cluster_codes=cluster_codes)
    new("daFrame", 
        SummarizedExperiment(
            assays=es,
            rowData=row_data,
            colData=col_data,
            metadata=metadata))
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
