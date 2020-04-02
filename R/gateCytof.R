#' @rdname gateCytof
#' @title Gating on CyTOF data
#' 
#' @description Wrapper around \code{openCyto} gating methods for 
#' \code{SingleCellExperiment}s, including group-specific gating parameters.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param chs character string of length 2 specifying the channels to gate on. 
#'   Should be one of \code{rownames(x)} or a \code{[int_]colData(x)} column,
#'   and correspond to a numeric variable.
#' @param type character string specifying the gate type:
#'   "rect" for rectangular, "elip" for eliptical, or "live" for polygonal.
#' @param xy for "elip" gates, a numeric vectors of length 2
#'   giving the x- and y- coordinates for the ellipse's center;
#'   for "rect" gates, a list of rectangular gate boundaries
#'   formatted as `list(c(xmin, ymin), c(xmax, ymax))`.
#' @param q numeric in (0,1) giving the quantile(s) 
#'   for gate of type "elip" and "live".
#' @param k integer number of cluster centers for eliptical gates.
#' @param i,s numeric specifying line intercept & slope for "live" gates.
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
#' @return a \code{SingleCellExperiment} with gating information stored in
#' \code{int_metadata(.)$gates$<gate_id>} as a list containing:
#' \describe{
#' \item{\code{chs}}{channels that were gated on (length 2 character vector)}
#' \item{\code{type}}{gating type (one of "rect", "elip", "live")}
#' \item{\code{pars}}{list of applied gating parameters
#' (e.g., group-specific quantiles for "elip"/"live" gates)}
#' \item{\code{gate_on}}{gate ID of an upstream gate determining a cell subset}
#' \item{\code{group_by}}{variable by which cells were grouped by}
#' \item{\code{data}}{\code{data.frame} storing the gate's metadata 
#' (e.g., x- and y-range for rectangular gates; 
#' gate boundaries for "elip" and "live" gates)}
#' }
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
#'   group_by = "file_name", gate_id = "cells",
#'   type = "elip", q = 0.95, xy = c(4, 4))
#'   
#' # view number & fraction of selected cells
#' table(sce$cells)
#' mean(sce$cells)
#' 
#' # view gating metadata
#' library(SingleCellExperiment)
#' gi <- int_metadata(sce)$gates
#' head(gi$cells$data) # data.frame of gate boundaries
#' gi$cells$pars       # list of applied parameters
#' 
#' # visualize gate on scatter plot
#' plotScatter(sce, gate_id = "cells")
#' 
#' @importFrom methods is
#' @importFrom S4Vectors metadata
#' @importFrom cluster silhouette
#' @importFrom dplyr bind_rows
#' @importFrom flowCore flowFrame flowSet 
#'   filterList ellipsoidGate rectangleGate
#' @importFrom flowWorkspace GatingSet lapply
#'   gh_pop_get_indices gs_pop_get_gate
#' @importFrom SingleCellExperiment int_colData int_metadata int_metadata<-
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom stats kmeans
#' @importFrom openCyto gs_add_gating_method gate_flowClust_2d
#' @export

gateCytof <- function(x, 
    chs, type = c("rect", "elip", "live"), 
    xy = c(1, 1), q = 0.99, k = NULL, i = 1, s = 0.5,
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
        is.numeric(q), q > 0, q < 1, is.numeric(i), is.numeric(s),
        is.null(k) || is.numeric(k) && (length(k) == 1 & k == round(k)),
        is.null(gate_id) || is.character(gate_id) & length(gate_id) == 1,
        is.null(group_by) || is.character(group_by)
        & length(group_by == 1) & group_by %in% vars,
        is.null(gate_on) || is.character(gate_on) & length(gate_on) == 1 
        & gate_on %in% names(int_metadata(x)$gates),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        is.logical(overwrite), length(overwrite) == 1)
    for (u in c("q", "i", "s"))
        assign(u, .get_gate_pars(get(u), x, group_by, u))
    
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
        dfs <- split(df, df[[group_by]])
    } else dfs <- list(root = df)
    # select numeric variables only
    dfs <- lapply(dfs, function(u) 
        u[, vapply(as.list(u), is.numeric, logical(1))])
    
    # 'openCyto' gating --------------------------------------------------------
    if (guess_k <- type == "elip" && is.null(k))
        message("'k' required for eliptical gates ",
            "but unspecified;\nusing Silhouette width ", 
            "to estimate number of cluster centers.")
    names(ids) <- ids <- names(dfs)
    res <- lapply(ids, function(id) {
        # check validity of 'xy' argument
        if (type != "live") 
            stopifnot(is(xy, ifelse(type == "rect", "list", "numeric")))
        if (type == "rect") 
            stopifnot(length(unlist(xy)) == 4)
        # construct'flowSet' 
        es <- as.matrix(dfs[[id]])
        fs <- flowSet(flowFrame(es))
        if (type == "live") {
            # register live-gate
            .live_gate(x, q[id], i[id], s[id])
        } else if (guess_k) {
            # estimate number of clusters using silhouette info
            # run k-means clustering for k = 1, ..., 10
            u <- dfs[[id]][, chs]
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
                k, q[id], paste(xy, collapse = ",")),
            live = NA)
        # construct 'GatingSet' & apply gate
        suppressMessages({
            gs <- GatingSet(fs)
            gs_add_gating_method(
                gs, alias = gate_id, 
                gating_method, gating_args,
                pop = "+",  parent = "root",
                dims = paste(chs, collapse = ","))
        })
        dat <- gs_pop_get_gate(gs, gate_id); names(dat) <- id
        list(dat = dat, ids = gh_pop_get_indices(gs, gate_id))
    })
    # store gate assignments in cell metadata
    idx <- unlist(map(dfs, "cell_id"))
    x[[gate_id]] <- FALSE
    x[[gate_id]][idx] <- unlist(map(res, "ids"))
    # get relevant gating parameters
    pars <- switch(type,
        rect = xy,
        live = list(q = q, xy = xy),
        elip = list(q = q, xy = xy, k = k))
    # store gating information internally
    gi <- lapply(ids, function(id) 
        .get_gate(res[[id]]$dat, type, group_by, q = q[id]))
    gi <- list(
        chs = chs, 
        type = type, 
        pars = pars,
        gate_on = gate_on, 
        group_by = group_by,
        data = bind_rows(gi))
    int_metadata(x)$gates[[gate_id]] <- gi
    # return SCE
    return(x)
}