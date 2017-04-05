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
            return(message("'exprs' and 'normed_bcs' should contain ",
                "the same number of events."))
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
