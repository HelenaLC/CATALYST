#' @rdname guessPanel
#' @title Guess parameter panel
#' 
#' @description Helper function to parse information from the 
#' \code{parameters} slot of a \code{flowFrame}/\code{flowSet}.
#'
#' @param x 
#'   a \code{\link[flowCore]{flowFrame}} or \code{\link[flowCore]{flowSet}}.
#' @param index 
#'   numeric. If \code{x} is a \code{flowSet} object, 
#'   this index specifies which \code{flowFrame} to extract.
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
#' @author Mark D Robinson, 
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @examples
#' # examplary data with Time, DNA, BC channels, etc.
#' data(raw_data)
#' guessPanel(raw_data)
#' 
#' @importFrom flowCore flowFrame flowSet
#' @importFrom methods is
# ------------------------------------------------------------------------------

setMethod(f="guessPanel", 
    signature=signature(x="flowFrame"), 
    definition=function(x) {
        # assume 'name' and 'desc' columns exist in @parameters
        param_df <- as(parameters(x), "data.frame")[, c("name", "desc")]
        colnames(param_df)[1] <- "fcs_colname"
        rownames(param_df) <- NULL
        
        # make some guesses of how to parse / what columns to use
        if ( any(grepl("_", param_df$desc)) ) {
            ss <- strsplit(param_df$desc,"_")
            param_df$desc.1 <- sapply(ss, .subset, 1)
            param_df$antigen <- sapply(ss, .subset, 2)
            param_df$use_channel <- !is.na(param_df$antigen)
            dont_use <- sapply(c("^BC", "dead", "DNA"), grep, param_df$antigen)
            param_df$use_channel[unlist(dont_use)] <- 0
        } else {
            param_df$antigen <- param_df$desc
            param_df$use_channel <- 1
        }
        param_df
    }
)

#' @rdname guessPanel
setMethod(f="guessPanel", 
    signature=signature(x="flowSet"), 
    definition=function(x, index=1) {  
        guessPanel(x[[index]])
    }
)
    