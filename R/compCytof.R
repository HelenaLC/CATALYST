# ==============================================================================
# Compensate CyTOF experiment
# ------------------------------------------------------------------------------

#' @rdname compCytof
#' @title Compensate CyTOF experiment
#' 
#' @description 
#' Compensates a mass spectrometry based experiment using a provided spillover
#' matrix, assuming a linear spillover in the experiment.
#'
#' @param x       
#' a \code{\link{flowFrame}} OR a character string specifying 
#' the location of FCS files that should be compensates.
#' @param y 
#' a spillover matrix.
#' @param out_path
#' a character string. If specified, compensated FCS files will be generated 
#' in this location. If \code{x} is a character string, file names will be 
#' inherited from uncompensated FCS files and given extension "_comped".
#' Defaults to NULL. 
#' @param method
#' one of "flow" or "nnls".
#' 
#' @details
#' If the spillover matrix (SM) does not contain the same set of columns as 
#' the input experiment, it will be adapted according to the following rules:
#' \enumerate{
#' \item{columns present in the SM but not in the input data 
#' will be removed from it}
#' \item{non-metal columns present in the input but not in the SM 
#' will be added such that they do neither receive nor cause spill}
#' \item{metal columns that have the same mass as a channel present in the SM 
#' will receive (but not emit) spillover according to that channel}
#' \item{if an added channel could potentially receive spillover (as it has 
#' +/-1M or +16M of, or is of the same metal type as another channel measured), 
#' a warning will be issued as there could be spillover interactions that have
#' been missed and may lead to faulty compensation}
#' }
#' 
#' @return 
#' Compensates the input \code{\link{flowFrame}} or, 
#' if \code{x} is a character string, all FCS files in the specified location. 
#' If \code{out_path=NULL} (the default), returns a \code{\link{flowFrame}} 
#' containing the compensated data. Otherwise, compensated data will be written 
#' to the specified location as FCS 3.0 standard files. 
#' 
#' @examples
#' # get single-stained control samples
#' data(ss_exp)
#' 
#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#' 
#' # debarcode
#' re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' spillMat <- computeSpillmat(x = re)
#' compCytof(x = ss_exp, y = spillMat)
#'
#' @author 
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' and Vito Zanotelli \email{vito.zanotelli@uzh.ch}
#' @importFrom flowCore flowFrame colnames exprs compensate
#' @importFrom nnls nnls
#' @importFrom stats setNames
#' @export
# ------------------------------------------------------------------------------

setMethod(f="compCytof",
    signature=signature(x="flowFrame", y="matrix"),
    definition=function(x, y, out_path=NULL, method="flow") {
        
        sm <- adaptSpillmat(y, flowCore::colnames(x))
        if (method == "flow") { 
            ff_comped <- flowCore::compensate(x, sm)
        } else if (method == "nnls") {
            es_comped <- t(apply(flowCore::exprs(x), 1, 
                function(row) nnls(t(sm), row)$x))
            ff_comped <- x
            colnames(es_comped) <- colnames(flowCore::exprs(x))
            rownames(es_comped) <- rownames(flowCore::exprs(x))
            flowCore::exprs(ff_comped) <- es_comped
        } else {
            stop("'method' should be one of \"flow\" or \"nnls\".")
        }
        
        if (!is.null(out_path)) {
            outNm <- file.path(out_path, 
                paste0(flowCore::identifier(x), "_comped.fcs"))
            suppressWarnings(flowCore::write.FCS(ff_comped, outNm))
        } else {
            ff_comped
        }
    })

# ------------------------------------------------------------------------------
#' @rdname compCytof
setMethod(f="compCytof",
    signature=signature(x="character", y="matrix"),
    definition=function(x, y, out_path=NULL, method="flow") {
        
        if (!file.exists(x))
            stop("x is neither a flowFrame nor a valid file/folder path.")
        fcs <- list.files(x, ".fcs", ignore.case=TRUE, full.names=TRUE)
        if (length(fcs) == 0)
            stop("No FCS files found in specified location.")
        ffs <- lapply(fcs, flowCore::read.FCS)
        
        if (is.null(out_path)) {
            lapply(ffs, function(i) compCytof(i, y, out_path, method))
        } else {
            out_nms <- gsub(x, out_path, 
                gsub(".fcs", "_comped.fcs", ignore.case=TRUE, fcs))
            for (i in seq_along(ffs)) {
                comped <- compCytof(ffs[[i]], y, out_path, method)
                suppressWarnings(flowCore::write.FCS(comped, out_nms[i]))
            }
        }
    })

# ------------------------------------------------------------------------------
#' @rdname compCytof
setMethod(f="compCytof",
    signature=signature(x="ANY", y="data.frame"),
    definition=function(x, y, out_path=NULL, method="flow") {
        compCytof(x, as.matrix(y), out_path, method)
    })