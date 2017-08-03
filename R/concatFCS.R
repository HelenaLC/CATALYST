# ==============================================================================
# concatinates all FCS files in the specified location
# ------------------------------------------------------------------------------

#' @rdname concatFCS
#' @title FCS file concatination
#' 
#' @description 
#' Concatinates all input data.
#'
#' @param x 
#' can be either a \code{\link{flowSet}}, a list of \code{\link{flowFrame}}s, 
#' a character specifying the location of the FCS files to be concatinated, 
#' or a vector of FCS file names.
#' @param out_path
#' an optional character string. If specified, an FCS file 
#' of the concatinated data will be written to this location. 
#' If NULL (default), a \code{\link{flowFrame}} will be returned.
#' @param by_time logical. 
#' Specifies whether files should be ordered by time of acquisition.
#' @param file_num logical. 
#' Specifies whether a file number column should be added.
#'
#' @return 
#' a \code{\link{flowFrame}} containing measurement intensities 
#' of all input data or a character of the FCS file name.
#' 
#' @examples
#' data(raw_data)
#' concatFCS(raw_data)
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom flowCore colnames exprs flowFrame flowSet fsApply isFCSfile
#' @importFrom methods as
#' @export

# ------------------------------------------------------------------------------

setMethod(f="concatFCS",
    signature=signature(x="flowSet"),
    definition=function(x, out_path=NULL, by_time=TRUE, file_num=FALSE) {
        
        n <- length(x)
        if (by_time) {
            # order by time
            bts <- keyword(x, "$BTIM")
            o <- order(bts)
            x <- x[o]
        } else {
            o <- seq_len(n)
        }
        nPars <- ncol(x[[1]])
        nEvents <- as.numeric(keyword(x, "$TOT"))
        
        # concatenate
        es <- fsApply(x, exprs)
        timeCol <- grep("time", colnames(x), ignore.case=TRUE)
        start <- c(1, cumsum(nEvents)+1)
        end <- start[-1]-1
        for (i in seq_along(x)[-1]) {
            inds <- start[i]:end[i]
            es[inds, timeCol] <- es[inds, timeCol]+es[end[i-1], timeCol]
        }
        
        # get descriptions
        d <- description(x[[1]])
        d$`$TOT` <- sum(nEvents)
        d$`GUID` <- gsub("_\\d+.fcs", "_concat.fcs", 
            d$ORIGINALGUID, ignore.case=TRUE)
        d <- d[!names(d) %in% c("$FIL", "FILENAME", "ORIGINALGUID")]
        d$`FILE LIST` <- paste(paste0(seq_len(n), 
            "-", keyword(x, "ORIGINALGUID")), collapse=",")
        d$`$BTIM` <- d$`$BTIM`
        d$`$ETIM` <- description(x[[n]])$`$ETIM`
        d$`$COM` <- "FCS files concatenated by CATALYST"
        d[paste0("$P", seq_len(nPars), "R")] <- df$range
        d$transformation <- "applied"
        
        # add file number column
        if (file_num) {
            nPars <- nPars+1
            es <- as.matrix(cbind(es, "FileNum"=
                    rep(seq_len(n)[o], nEvents[o])))
            d[paste0("$P", nPars, "B")] <- 32
            d[paste0("$P", nPars, "E")] <- 0
            d[paste0("$P", nPars, "N")] <- "FileNum"
            d[paste0("$P", nPars, "R")] <- n-1
            d[paste0("$P", nPars, "S")] <- "File Number"
        }
        
        # get parameters
        pars <- colnames(es)
        PnS <- d[paste0("$P", seq_len(nPars), "S")]
        PnS[sapply(PnS, is.null)] <- NA
        
        mins <- matrixStats::colMins(es)
        maxs <- matrixStats::colMaxs(es)
        df <- data.frame(list(
            name=colnames(es), 
            desc=unlist(PnS), 
            range=abs(maxs)-abs(mins),
            minRange=mins,
            maxRange=maxs))
        p <- Biobase::AnnotatedDataFrame(df)
        rownames(p) <- paste0("$P", seq_len(nPars))

        # construct new flowFrame
        ff <- flowFrame(exprs=es, parameters=p, description=d)

        if (is.null(out_path)) return(ff)
        suppressWarnings(flowCore::write.FCS(ff, 
            file.path(out_path, description(ff)$GUID)))
    })

# ------------------------------------------------------------------------------

#' @rdname concatFCS
setMethod(f="concatFCS",
    signature=signature(x="character"),
    definition=function(x, out_path=NULL, by_time=TRUE, file_num=FALSE) {
        if (length(x) == 1) {
            fcs <- list.files(x, ".fcs", full.names=TRUE, ignore.case=TRUE)
        } else {
            fcs <- x
        } 
        n <- length(fcs)
        if (sum(flowCore::isFCSfile(fcs)) != n) 
            stop("Not all files in ", x, " are valid FCS files.")   
        if (n < 2) {
            if (n == 0) 
                stop("No FCS files have been found in \"", x, "\"")
            if (n == 1) 
                stop("Only a single FCS file has been found in \"", x,"\"") 
        }
        ffs <- lapply(fcs, flowCore::read.FCS, 
            transformation=FALSE, truncate_max_range=FALSE)
        concatFCS(ffs, out_path, by_time, file_num)
    })

# ------------------------------------------------------------------------------

#' @rdname concatFCS
setMethod(f="concatFCS",
    signature=signature(x="list"),
    definition=function(x, out_path=NULL, by_time=TRUE, file_num=FALSE) {
        if (length(x) == 1) 
            stop("Only a single flowFrame has been provided.")
        concatFCS(as(x, "flowSet"), out_path, by_time, file_num)
    })