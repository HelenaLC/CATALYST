#' @rdname guessPanel
#' @title Guess parameter panel
#' 
#' @description Helper function to parse information from the 
#'   \code{parameters} slot of a \code{flowFrame}/\code{flowSet}.
#'
#' @param x a \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param sep character string specifying how channel descriptions 
#'   should be parsed. E.g., if \code{pData(x)$desc} contains both 
#'   channel and antigens formatted as, \code{155Gd_CD73}, 
#'   descriptions will be split according to \code{sep} and 
#'   everything after the first \code{sep} will be used as 
#'   the antigen name (here, CD73).
#'
#' @return a \code{data.frame} with the following columns:
#' \itemize{
#' \item{\code{name}: the parameter name 
#' as extracted from the input \code{flowFrame},}
#' \item{\code{desc}: the parameter description 
#' as extracted from the input \code{flowFrame},}
#' \item{\code{antigen}: the targeted protein markers, and}
#' \item{\code{use_channel}: logical. If TRUE, the channel 
#' is expected to contain a marker and will be kept.}}
#' 
#' @author Mark D Robinson &
#' Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @examples
#' # examplary data with Time, DNA, BC channels, etc.
#' data(raw_data)
#' guessPanel(raw_data[[1]])
#' 
#' @importFrom dplyr mutate_if rename
#' @importFrom flowCore parameters
#' @importFrom methods as is
#' @importFrom purrr map
#' @export

guessPanel <- function(x, sep = "_") {
    # check validity of input arguments
    stopifnot(is(x, "flowFrame"),
        is.character(sep), length(sep) == 1)
    ps <- as(parameters(x), "data.frame")
    ps <- mutate_if(ps, is.factor, as.character)
    stopifnot((k <- c("name", "desc")) %in% colnames(ps))
    ps <- ps[, k]

    # do this if descriptions contain channel masses
    ms1 <- .get_ms_from_chs(ps$name)
    ms2 <- .get_ms_from_chs(ps$desc)
    check <- map(lapply(ms1[!is.na(ms1)], grep, x = ms2), 1)
    if (length(unlist(check)) == sum(!is.na(ms1))) {
        # make some guesses of how to parse / what columns to use
        if (any(grepl(sep, ps$desc))) {
            # split on 1st occurange of 'sep'
            ss <- strsplit(ps$desc, sep)
            ss <- lapply(ss, function(u) if (length(u) > 2) 
                c(u[1], paste(u[-1], collapse = sep)) else u)
            ps <- rename(ps, desc0 = "desc")
            ps$desc <- map(ss, 1)
            ps$antigen <- map(ss, function(u)
                ifelse(length(u) == 1, u[1], paste(u[-1], collapse = sep)))
            ps$use_channel <- !is.na(ps$antigen)
            ps$use_channel[paste(ps$desc) == paste(ps$antigen)] <- FALSE
        } else {
            ps$antigen <- ps$desc
            ps$use_channel <- TRUE
        }
        ex <- c("^(BC|MCB)[0-9]", "dead", "dna", "bead", "vol")
        dont_use <- lapply(ex, grep, ps$antigen, ignore.case = TRUE)
        ps$use_channel[unlist(dont_use)] <- FALSE
    } else {
        # use available descriptions otherwise
        ps$antigen <- ps$desc
        ps$use_channel <- !is.na(ps$antigen)
    }
    rownames(ps) <- NULL
    colnames(ps)[1] <- "fcs_colname"
    return(ps)
}
