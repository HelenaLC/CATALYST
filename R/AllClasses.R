# ================================================================================
# Debarcoding frame class
# --------------------------------------------------------------------------------

#' @rdname dbFrame-class
#' @name dbFrame-class
#' 
#' @title Debarcoding frame class
#' @description Class representing the data returned by and used throughout debarcoding.
#' 
#' @param object a \code{dbFrame}.
#' 
#' @details 
#' \code{show(dbFrame)} will display \itemize{
#' \item the dimensionality of the measurement data and number of barcodes used
#' \item current assignments in order of decreasing population size
#' \item separation cutoff estimates
#' \item the average and per-population yield achieven upon debarcoding}
#' 
#' @slot exprs  
#' a matrix containing raw intensities of the input flowFrame.
#' @slot bc_key 
#' binary barcoding scheme with numeric masses as column names OR 
#' a numeric vector specifying indices of columns corresponding to barcode channels.
#' @slot bc_ids
#' numeric vector of barcode IDs. If a barcoding scheme is supplied, the respective 
#' binary code's row index, else, the column index of the respective barcode channel.
#' @slot deltas 
#' numeric vector of separations between positive and negative barcode populations 
#' computed from normalized barcode intensities.
#' @slot normed_bcs 
#' matrix containing normalized barcode intensities.
#' @slot sep_cutoffs
#' numeric vector of distance separation cutoffs between positive and negative 
#' barcode populations above which events will be unassigned.
#' @slot mhl_cutoff
#' non-negative numeric value specifying the Mahalanobis distance cutoffs
#' below which events will be unassigned.
#' @slot counts
#' matrix of dimension (# barcodes)x(101) where each row contains
#' the number of events within a barcode for which positive and negative populations 
#' are separated by a distance between in [0,0.01), ..., [0.99,1], respectively.
#' @slot yields
#' a matrix of dimension (# barcodes)x(101) where each row contains
#' the percentage of events within a barcode that will be obtained after applying
#' a separation cutoff of 0, 0.01, ..., 1, respectively.
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @importFrom methods new
#' @export

# --------------------------------------------------------------------------------

dbFrame <- setClass(Class="dbFrame", 
    package="CATALYST",
    slots=c(exprs="matrix",
        bc_key="data.frame",
        bc_ids="numeric",
        deltas="numeric",
        normed_bcs ="matrix",
        sep_cutoffs="numeric",
        mhl_cutoff="numeric",
        counts="matrix",
        yields="matrix"))

# --------------------------------------------------------------------------------

setValidity(Class="dbFrame", 
    method=function(object){
        if(!any(object@bc_ids %in% c(0, as.numeric(rownames(object@bc_key)))))
            return(cat("Invalid 'bc_ids': Barcode IDs must be numeric",
                "\nand occur in the row names of 'bc_key'."))
        
        nms <- colnames(object@exprs)
        ms <- as.numeric(regmatches(nms, gregexpr("[0-9]+", nms)))
        if(!all(colnames(object@bc_key) %in% ms))
            return(cat("Invalid 'bc_key': Column names must be numeric",
                "\nand coherent with masses extracted from 'exprs'."))
        
        return(TRUE)
    })










