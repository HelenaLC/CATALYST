# ==============================================================================
# concatenation expression matrices of a flowSet
# ------------------------------------------------------------------------------

.concat_fs <- function(fs, ns) {
    es <- fsApply(fs, exprs)
    timeCol <- grep("time", colnames(fs), ignore.case=TRUE)
    start <- c(1, cumsum(ns)+1)
    end <- start[-1]-1
    for (i in seq_along(fs)[-1]) {
        inds <- start[i]:end[i]
        es[inds, timeCol] <- es[inds, timeCol]+es[end[i-1], timeCol]
    }
    return(es)
}
