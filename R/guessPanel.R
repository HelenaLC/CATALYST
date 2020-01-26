#' @rdname guessPanel
#' @title Guess parameter panel
#' 
#' @description Helper function to parse information from the 
#'   \code{parameters} slot of a \code{flowFrame}/\code{flowSet}.
#'
#' @param x a \code{\link[flowCore]{flowSet}}.
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
#' @export

guessPanel <- function(x) {
    stopifnot(is(x, "flowFrame"))
    ps <- as(parameters(x), "data.frame")
    ps <- mutate_if(ps, is.factor, as.character)
    stopifnot((k <- c("name", "desc")) %in% colnames(ps))
    ps <- ps[, k]
    rownames(ps) <- NULL
    colnames(ps)[1] <- "fcs_colname"

    # make some guesses of how to parse / what columns to use
    if (any(grepl("_", ps$desc))) {
        ss <- strsplit(ps$desc, "_")
        ps <- rename(ps, desc0 = "desc")
        ps$desc <- vapply(ss, .subset, i = 1, character(1))
        ps$antigen <- vapply(ss, .subset, i = 2, character(1))
        ps$use_channel <- !is.na(ps$antigen)
        dont_use <- lapply(c("^BC", "dead", "DNA", "B"), grep, ps$antigen)
        ps$use_channel[unlist(dont_use)] <- FALSE
    } else {
        ps$antigen <- ps$desc
        ps$use_channel <- TRUE
    }
    ps
}