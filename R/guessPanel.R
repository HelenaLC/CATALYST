#' Helper function to parse information from the @parameters slot of a flowFrame/flowSet object
#'
#' @param ff a flowSet or flowFrame object
#' @param index if ff is a flowSet object, the index is used to extract a flowFrame
#'
#' @return
#' @export
#'
#' @examples
guessPanel <- function(ff, index = 1) {
  isFS <- is(ff, "flowSet")
  stopifnot( is(ff, "flowFrame") | isFS )
  if(isFS)
    ff <- ff[[1]]

  # assume 'name' and 'desc' columns exist in @parameters
  param_df <- as(ff@parameters,"data.frame")[,c("name","desc")]
  colnames(param_df)[index] <- "fcs_colname"
  rownames(param_df) <- NULL

  # make some guesses of how to parse / what columns to use
  if( any(grepl("_", param_df$desc)) ) {
      ss <- strsplit(param_df$desc,"_")
      param_df$desc.1 <- sapply(ss, .subset, 1)
      param_df$antigen <- sapply(ss, .subset, 2)
      param_df$use_channel <- !is.na(param_df$antigen)+0
      param_df$use_channel[grep("^BC",param_df$antigen)] <- 0
      param_df$use_channel[grep("dead",param_df$antigen)] <- 0
      param_df$use_channel[grep("DNA",param_df$antigen)] <- 0
  } else {
      param_df$antigen <- param_df$desc
      param_df$use_channel <- 1
  }
  param_df
}

