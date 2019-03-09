#' @rdname daFrame-methods
#' @param i,j indices specifying elements to extract.
#' @importFrom dplyr %>% mutate_if
#' @importFrom methods as
#' @importFrom SingleCellExperiment int_elementMetadata SingleCellExperiment
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom S4Vectors metadata SimpleList
#' @export

setMethod("[", 
    c("daFrame", "ANY", "ANY"),
    function(x, i, j) {
        if (missing(i)) i <- TRUE else if (is.logical(i)) 
            stopifnot(length(i) %in% c(1, nrow(x)))
        if (missing(j)) j <- TRUE else if (is.logical(j)) 
            stopifnot(length(j) %in% c(1, ncol(x)))
        # store original reducedDims
        dr <- reducedDims(x)
        is <- seq_len(nrow(x))[i]
        int_em <- int_elementMetadata(x)
        reducedDims(x) <- NULL
        # subset assays, colData, and rowData
        y <- as(x, "SingleCellExperiment")
        y <- y[i, j]
        # drop missing rowData levels
        rowData(y) <- rowData(y) %>% 
            data.frame %>% 
            mutate_if(is.factor, droplevels)
        # update metadata
        x <- as(y, "daFrame")
        ei <- metadata(x)$experiment_info
        sids <- sample_ids(x)
        samples_keep <- ei$sample_id %in% levels(sids)
        ei <- ei[samples_keep, ] %>% 
            mutate_if(is.factor, droplevels)
        m <- match(ei$sample_id, levels(sids))
        ei$n_cells <- as.numeric(table(sids)[m])
        metadata(x)$experiment_info <- ei
        # subset reducedDims
        drs <- names(dr)
        names(drs) <- drs
        dr <- lapply(drs, function(u) {
            foo <- matrix(nrow = 1, ncol = nrow(dr[[u]]))
            foo <- SingleCellExperiment(foo, reducedDims = dr[u])
            k <- grep(u, names(int_em))
            idx <- which(int_em[[k]])
            is <- intersect(is, idx)
            m <- match(idx, is, nomatch = 0)
            out <- reducedDim(foo[, m])
            attr(out, "percentVar") <- attr(dr[[u]], "percentVar")
            return(out)
        }) %>% SimpleList
        drop <- vapply(dr, nrow, numeric(1)) == 0
        if (any(drop))
            message("Dropping DR(s) out of bounds: ", 
                paste(dQuote(names(dr)[drop]), collapse = ", "))
        x@reducedDims <- dr[!drop]
        return(x)
    }
)
