#' @rdname gateCytof
#' @title Gating on CyTOF data
#' 
#' @description 
#' Wrapper around dimension reduction methods available 
#' through \code{scater}, with optional cell subsampling.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param chs character string of length 2 specifying the channels to gate on. 
#'   Should be one of \code{rownames(x)} or a \code{[int_]colData(x)} column,
#'   and correspond to a numeric variable.
#' @param q numeric in (0,1) giving the quantile for eliptical and live-gating.
#' @param k integer number of cluster centers for eliptical gating.
#' @param bs numeric of length 2 specifying the intercept & slope
#'   of the line used for live-gating.
#' @param gate_id character string giving a unique identifier for the gate.
#' @param group_by character string specifying the grouping variable
#'   in case gates should be applied to specified cell subsets (e.g., 
#'   indenpendently for each sample). Should be a \code{[int_]colData(x)} 
#'   column, and correspond to a logical or factor variable.
#' @param gate_on character string specifying a previously applied gate
#'   to gate on, i.e., from which to subset cells prior to gating.
#'   Should be one of \code{names(int_metadata(x)$gates)}.
#' @param assay character string specifying which assay data to use.
#'   Should be one of \code{assayNames(x)}.
#' @param overwrite logical. If there already exists a gate 
#'   with identifier \code{gate_id}, should it be overwritten?
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @return a \code{SingleCellExperiment}.
#' 
#' @examples
#' data(raw_data)
#' sce <- fcs2sce(raw_data)
#' 
#' # GATING FOR CELLS
#' 
#' # specify DNA channels
#' dna_chs <- c("Ir191Di", "Ir193Di")
#' 
#' # apply eliptical gate
#' sce <- gateCytof(sce, dna_chs, 
#'   type = "elip", q = 0.95, xy = c(4, 4), 
#'   gate_id = "cells")
#'   
#' # view number & fraction of selected cells
#' table(sce$cells)
#' mean(sce$cells)
#' 
#' # visualize gate on scatter plot
#' plotScatter(sce, gate_id = "cells")
#' 
#' @importFrom methods is
#' @importFrom dplyr select_if
#' @importFrom S4Vectors metadata
#' @importFrom cluster silhouette
#' @importFrom flowCore flowFrame flowSet 
#'   filterList ellipsoidGate rectangleGate
#' @importFrom flowWorkspace GatingSet lapply
#'   gh_pop_get_indices gs_pop_get_gate
#' @importFrom SingleCellExperiment int_colData int_metadata int_metadata<-
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom stats kmeans
#' @importFrom openCyto gs_add_gating_method gate_flowClust_2d
#' @export

gateCytof <- function(x, chs, 
    xy = NULL, q = 0.99, k = NULL, bs = c(1, 0.5),
    type = c("rect", "quad", "elip", "live"), 
    gate_id = NULL, group_by = NULL, gate_on = NULL,
    assay = "exprs", overwrite = FALSE) {
    # check validity of input arguments
    type <- match.arg(type)
    stopifnot(is(x, "SingleCellExperiment"))
    # extract valid variables (channels & any cell metadata)
    vars <- c(rownames(x), names(colData(x)), names(int_colData(x)))
    stopifnot(
        is.null(xy) || is.numeric(unlist(xy)) & length(xy) == 2,
        is.character(chs), length(chs) == 2, chs %in% vars,
        is.numeric(q), length(q) == 1, q > 0, q < 1,
        is.null(k) || is.numeric(k) && (length(k) == 1 & k == round(k)),
        is.numeric(bs), length(bs) == 2,
        is.null(gate_id) || is.character(gate_id) & length(gate_id) == 1,
        is.null(group_by) || is.character(group_by)
        & length(group_by == 1) & group_by %in% vars,
        is.null(gate_on) || is.character(gate_on) & length(gate_on) == 1 
        & gate_on %in% names(int_metadata(x)$gates),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        is.logical(overwrite), length(overwrite) == 1)
    # construct unique ID if 'gate_id' is not specified
    if (is.null(gate_id)) {
        n_gs <- length(int_metadata(x)$gates)
        gate_id <- paste0("gate", n_gs + 1)
     # when 'overwrite = FALSE', check that gate ID doesn't exist already 
    } else if (!overwrite && gate_id %in% c(
        names(colData(x)), names(int_metadata(x)$gates))) 
        stop("Gate with ID ", dQuote(gate_id), " alread exists;\n ", 
            " specify another 'gate_id' or 'overwrite = TRUE'.")
    # construct data.frame containing specified 
    # assay data & all numeric cell metadata
    cd <- cbind(colData(x), int_colData(x))
    cd <- cd[, vapply(cd@listData, class, character(1)) != "numeric"]
    df <- data.frame(t(as.matrix(assay(x, assay))), cd,
        check.names = FALSE, stringsAsFactors = FALSE,
        # add cell ID column to track cells
        cell_id = seq_len(ncol(x))) 
    # subset cells if 'gate_on' is specified
    if (!is.null(gate_on))
        df <- df[x[[gate_on]], ]
    # split cells if 'group_by' is specified
    if (!is.null(group_by)) {
        stopifnot(!is.numeric(df[[group_by]]))
        df <- split(df, df[[group_by]])
    } else df <- list(df)
    # select numeric variables only
    df <- lapply(df, select_if, is.numeric)
    # 'openCyto' gating methods
    if (type != "quad") {
        if (type != "live") stopifnot(is(xy, 
            ifelse(type == "rect", "list", "numeric")))
        if (type == "rect") stopifnot(length(unlist(xy)) == 4)
        # construct'flowSet' 
        ffs <- lapply(df, function(u) 
            flowFrame(as.matrix(u)))
        fs <- flowSet(ffs)
        if (type == "live") {
            # register live-gate
            .live_gate(x, q, bs)
        } else if (type == "elip" && is.null(k)) {
            message("'k' required for eliptical gates ",
                "but unspecified;\nusing Silhouette width ", 
                "to estimate number of cluster centers.")
            # estimate number of clusters using silhouette info
            # run k-means clustering for k = 1, ..., 10
            u <- df[[1]][, chs]
            d <- dist(u)
            names(ks) <- ks <- seq_len(5)
            km_res <- lapply(ks, kmeans, x = u)
            # compute mean Silhouette score for each k
            scores <- vapply(ks[-1], function(k) {
                kids <- km_res[[k]]$cluster
                res <- silhouette(kids, d)
                mean(res[, "sil_width"])
            }, numeric(1))
            # use 'k = 1' if all bad, otherwise use best scoring k
            k <- ifelse(!any(scores > 0.8), 1, which.max(scores) + 1)
            message("=> Using 'k = ", k, "'.")
        }
        gating_method <- switch(type,
            rect = "boundary",
            elip = "flowClust.2d",
            live = "liveGate")
        gating_args <- switch(type,
            rect = sprintf(
                "min=c(%s),max=c(%s)", 
                paste(xy[[1]], collapse = ","),
                paste(xy[[2]], collapse = ",")),
            elip = sprintf(
                "K=%s,quantile=%s,target=c(%s)",
                k, q, paste(xy, collapse = ",")),
            live = NA)
        # construct 'GatingSet' & apply gate
        suppressMessages({
            gs <- GatingSet(fs)
            gs_add_gating_method(gs,
                alias = gate_id, pop = "+",  parent = "root", 
                dims = paste(chs, collapse = ","),
                gating_method, gating_args)
        })
        ids <- lapply(gs, gh_pop_get_indices, gate_id)
        gs <- gs_pop_get_gate(gs, gate_id)
    # quadrant-type gate
    } else {
        ids <- lapply(df, function(u) {
            # subset channels of interest
            u <- u[, chs]
            # test for right & upper half
            r <- u[, 1] > xy[1]
            u <- u[, 2] > xy[2]
            # test which quadrant ea. cell falls into
            ids <- list(r & u, !r & u, !r & !u, r & !u)
            ids <- do.call("cbind", ids)
            apply(ids, 1, which)
        })
    }
    # store gate assignments
    idx <- unlist(map(df, "cell_id"))
    x[[gate_id]] <- FALSE
    x[[gate_id]][idx] <- unlist(ids)
    # store gating information 
    gi <- list(type = type, chs = chs,
        group_by = group_by, gate_on = gate_on, 
        data = .get_gate(gs, type, group_by, q = q))
    int_metadata(x)$gates[[gate_id]] <- gi
    # return SCE
    return(x)
}