#' @rdname filter
#' @title Filter daFrame
#' 
#' @description 
#' Filters events from a \code{daFrame} using conditional statements.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param ... 
#'   conditional statements separated by comma.
#'   Left-hand side arguments must occur in the \code{colnames(rowData(x))};
#'   accepted operators are \code{==}, \code{!=} and \code{\%in\%}.
#' @param k 
#'   numeric or character string. 
#'   Specifies the clustering to extract populations from.
#'   Must be one of \code{names(cluster_codes(x))}.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @return a \code{daFrame}.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md, merging_table)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run clustering
#' re <- cluster(re)
#' 
#' # one condition only, remove a single sample
#' filter(re, condition == "Ref", sample_id != "Ref1")
#' 
#' # keep only a subset of clusters
#' filter(re, cluster_id %in% c(7, 8, 18), k = 20)
# ------------------------------------------------------------------------------

setMethod(f="filter", 
    signature=signature(x="daFrame"), 
    definition=function(x, ..., k = NULL) {
        
        n <- nrow(x)
        args <- substitute(deparse(...))
        args <- as.character(args)[-1]
        args <- sapply(args, function(x) gsub("\"", "\'", x, fixed = TRUE))
        
        # check validity of left- & right-hand arguments
        args_split <- sapply(args, strsplit, "(==).|(!=).|(%in%).")
        l <- sapply(args_split, "[[", 1)
        l <- sapply(l, function(x) gsub("(!)|[[:blank:]]", "", x))
        stopifnot(all(l %in% colnames(rowData(x))))
        
        r <- sapply(args_split, "[[", 2)
        md <- metadata(x)$experiment_info
        for (i in seq_along(args)[l != "cluster_id"])
            stopifnot(all(gsub("\'", "", r[[i]]) %in% md[, l[[i]]]))
        
        # filter events
        if (any(l != "cluster_id")) {
            for (i in colnames(rowData(x))) assign(i, rowData(x)[, i])
            inds <- sapply(args[l != "cluster_id"], 
                function(x) eval(parse(text=x)))
            inds <- apply(inds, 1, all)
            if (sum(inds) == 0) 
                stop("The applied filtering would remove all events.")
            x <- x[inds, ]
        }
        
        if ("cluster_id" %in% l) {
            if (is.null(k))
                stop("Please specify which clustering 'k' to use.")
            
            # check that cluster() has been run
            stopifnot("cluster_codes" %in% names(metadata(x)))
            stopifnot("cluster_id" %in% names(rowData(x)))
            
            # check validity of cluster IDs
            k <- check_validity_of_k(x, k)
            stopifnot(all(eval(parse(text=r[l == "cluster_id"])) %in% cluster_codes(x)[, k]))
            
            cluster_id <- cluster_codes(x)[cluster_ids(x), k]
            inds <- eval(parse(text = args[l == "cluster_id"]))
            if (sum(inds) == 0) 
                stop("The applied filtering would remove all events.")
            x <- x[inds, ]
        }
        
        # update factor levels
        for (i in seq_len(ncol(rowData(x))))
            rowData(x)[, i] <- factor(as.character(rowData(x)[, i]))
        
        # update metadata
        for (i in colnames(md)) assign(i, md[, i])
        if (any(l != "cluster_id")) {
            inds <- sapply(args[l != "cluster_id"], 
                function(x) eval(parse(text=x)))
            inds <- apply(inds, 1, all)
            md <- md[inds, ]
            rownames(md) <- NULL
            metadata(x)$experiment_info <- md
        }
        metadata(x)$n_cells <- table(sample_ids(x))
        message(n - nrow(x), " / ", n, " events have been removed.")
        return(x)
    }
)
