#' @rdname daFrame-methods
#' @param i,j indices specifying elements to extract.
#' @importFrom dplyr %>% mutate_if
#' @importFrom methods as
#' @importFrom SingleCellExperiment int_elementMetadata SingleCellExperiment
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom S4Vectors metadata
#' @export
setMethod("[", 
    c("daFrame", "ANY", "ANY"),
    function(x, i, j) {
        
        if (missing(i) && missing(j))
            return(x)
        
        if (missing(i) || is.logical(i) && isTRUE(i)) {
            i <- seq_len(nrow(x))
        } else if (is.logical(i)) {
            if (length(i) != 1) 
                stopifnot(length(i) == nrow(x))
            i <- which(i)
        }
        
        if (missing(j) || is.logical(j) && isTRUE(j)) {
            j <- seq_len(ncol(x))
        } else if (is.logical(j)) {
            if (length(j) != 1)
                stopifnot(length(j) == ncol(x))
            j <- which(j)
        } else if (is.character(j)) {
            m <- match(j, colnames(x))
            nomatch <- is.na(m)
            if (any(nomatch)) {
                nomatch <- unique(j[nomatch])
                stop("Index out of bonds: ",
                    paste(nomatch, collapse = ", "))
            }
            j <- m
        }
    
        # subset assays
        as <- lapply(assays(x), "[", i, j, drop = FALSE)   
        # subset rowData
        rd <- rowData(x)[i, ] %>% data.frame %>% 
            mutate_if(is.factor, droplevels)  
        # subset colData
        cd <- colData(x)[j, ]
        # subset reducedDims
        drs <- reducedDimNames(x)
        names(drs) <- drs
        dr <- lapply(drs, function(dr) {
            idx <- int_elementMetadata(x)
            k <- grep(dr, names(idx))
            idx <- which(idx[[k]])
            i <- intersect(i, idx)
            m <- match(idx, i, nomatch = 0)
            y <- reducedDim(x, dr)
            pv <- attr(y, "percentVar")
            y <- y[m, , drop = FALSE]
            attr(y, "percentVar") <- pv
            return(y)
        })
        dr <- as(dr, "SimpleList")
        int_em <- int_elementMetadata(x)
        if (length(int_em) != 0) {
            int_em <- int_em[i, , drop = FALSE]
        } else {
            int_em <- matrix(nrow = length(i), ncol = 0)
            int_em <- DataFrame(int_em)
        }
        # return daFrame
        sce <- SingleCellExperiment(
            assays = as,
            rowData = rd,
            colData = cd,
            metadata = metadata(x))
        x <- as(sce, "daFrame")
        x@reducedDims <- dr
        x@int_elementMetadata <- int_em
        return(x)   
    }
)