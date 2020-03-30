#' @rdname plotScatter
#' @title Scatter plot
#' 
#' @description 
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param chs character string of length 2 specifying which channels to plot.
#'   Should be one of \code{rownames(x)} or \code{names([int_]colData(x))},
#'   and correspond to a numeric or logical variable.
#' @param assay character string specifying which assay data to use.
#'   Should be one of \code{assayNames(x)}.
#' @param gate_id a character vector specifying 
#'   a gate added via \code{\link{gateCytof}}.
#' @param show_gate logical. If \code{gate_id} is specified, specifies 
#'   whether the boundary of the specified gate should be visualized.
#'   Ignored if `chs` don't correspond to those used by the gate.
#' @param show_perc logical. If \code{gate_id} is specified, specifies 
#'   whether the percentage of gated cells should be included.
#' @param color_by specifies the color coding. If \code{gate_id} is specified, 
#'   \code{color_by = "selection"} will color cells in a binary fashion, 
#'   i.e., whether cells fall within the gate's boundary or not.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(raw_data)
#' sce <- fcs2sce(raw_data)
#' 
#' # apply file-specific eliptical gate
#' sce <- gateCytof(sce, 
#'   chs = c("Ir191Di", "Ir193Di"),
#'   gate_id = "cells",
#'   group_by = "file_name",
#'   type = "elip", 
#'   xy = c(4, 4))
#' 
#' # visualize gates on scatter plots
#' plotScatter(sce, gate_id = "cells")
#' 
#' @import ggplot2
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom matrixStats colMaxs
#' @importFrom SingleCellExperiment int_colData
#' @importFrom SummarizedExperiment assay assayNames colData
#' @export

plotScatter <- function(x, chs = NULL, gate_id = NULL, 
    show_gate = TRUE, show_perc = TRUE, bins = 100,
    assay = "exprs", color_by = c("density", "selection")) {
    # check validity of input arguments
    color_by <- match.arg(color_by)
    stopifnot(is(x, "SingleCellExperiment"))
    cd_vars <- c(names(colData(x)), names(int_colData(x)))
    stopifnot(is.null(chs) & !is.null(gate_id) 
        || is.character(chs) & chs %in% c(rownames(x), cd_vars),
        is.character(assay), length(assay) == 1, assay %in% assayNames(x),
        is.null(gate_id) || is.character(gate_id) & length(gate_id) == 1 
        & !is.null(gi <- int_metadata(x)$gates[[gate_id]]))
    # construct data.frame of specified assay data & all cell metadata
    df <- data.frame(check.names = FALSE, stringsAsFactors = FALSE,
        t(as.matrix(assay(x, assay))), colData(x), int_colData(x))
    cd_vars <- intersect(cd_vars, names(df))
    # if 'gate_id' is specified, 
    # default to channels used by gate
    if (!is.null(gate_id) & is.null(chs)) chs <- gi$chs
    # get axis limits prior to potential subsetting
    maxs <- colMaxs(as.matrix(df[, chs]))
    maxs <- vapply(maxs, function(u) ceiling(u*2)/2, numeric(1))
    lims <- list(c(-0.5, maxs[1]), c(-0.5, maxs[2]))
    #lims <- c(-0.5, ceiling(max(df[, chs])*2)/2)
    # subset if parent != all cells
    if (!is.null(gate_id) && !is.null(gi$gate_on))
        df <- df[df[[gi$gate_on]], ]
    # initialize facetting & (optionally) melt data.frame 
    if (length(chs) > 2) {
        df <- melt(df, id.vars = c(chs[1], cd_vars))
        df <- df[df$variable %in% chs, ]
        facet <- facet_wrap("variable")
        ylab <- ylab(NULL)
        chs[2] <- "value"
    } else {
        facet <- NULL
        ylab <- NULL
    }
    # initial aesthetics
    fill_var <- "sqrt(..ncount..)"
    fill_pal <- scale_fill_distiller(palette = "Spectral")
    gate_geom <- geom_perc <- NULL
    if (is.null(gate_id)) {
        if (color_by == "selection")
            message("No 'gate_id' specified. Coloring by density instead.")
    } else {
        if (color_by == "selection") {
            fill_var <- gate_id
            fill_pal <- scale_fill_manual(values = c("grey", "royalblue"))
        } 
        # update facetting scheme
        if (is.null(facet)) {
            if (!is.null(gi$group_by))
                facet <- facet_wrap(gi$group_by)
        } else if (!is.null(gi$group_by)) {
            facet <- facet_grid(
                rows = vars(!!sym(gi$group_by)), 
                cols = vars(!!sym("variable")))
            # drop empty groups from gate data
            d <- gi$data; g <- gi$group_by
            gi$data <- d[d[[g]] %in% df[[g]], ]
        }
        # (optionally) add fraction of gates cells
        if (show_perc) {
            if (!is.null(gi$group_by)) {
                cs <- split(x[[gate_id]], x[[gi$group_by]])
                ps <- vapply(cs, mean, numeric(1))
            } else {
                ps <- mean(x[[gate_id]])
            }
            ps <- sprintf("%1.2f%%", ps*100)
            ps <- data.frame(p = ps)
            if (!is.null(gi$group_by))
                ps[[gi$group_by]] <- names(cs)
            geom_perc <- geom_label(
                data = ps, inherit.aes = FALSE, 
                label.r = unit(0.1, "lines"), 
                label.padding = unit(0.2, "lines"),
                hjust = 0, vjust = 0, col = "red",
                aes_string("-0.3", "-0.3", label = "p"))
        }
        # (optionally) add geom for gate
        if (show_gate) {
            check1 <- is.null(facet)
            check2 <- !is(facet, "FacetGrid") 
            check3 <- length(facet$params$facets) != 2
            if (check1 || isTRUE(check2 && check3)) {
                gate_geom <- switch(gi$type, 
                    rect = geom_rect(
                        data = gi$data, aes_string(
                            xmin = "xmin", xmax = "xmax", 
                            ymin = "ymin", ymax = "ymax"),
                        inherit.aes = FALSE, col = "red", fill = NA),
                    elip = geom_path(
                        data = gi$data, aes_string("x", "y"),
                        inherit.aes = FALSE, col = "red"),
                    live = geom_path(
                        data = gi$data, aes_string("x", "y"),
                        inherit.aes = FALSE, col = "red"))
            } else {
                message("Specified 'chs' don't match with those used by",
                    " gate ", dQuote(gate_id), "; ignoring 'show_gate'.")
            }
        }
    } 
    xy <- sprintf("`%s`", chs)
    ggplot(df, aes_string(xy[1], xy[2], fill = fill_var)) + 
        geom_hex(na.rm = TRUE, bins = bins) + 
        fill_pal + guides(fill = FALSE) + 
        scale_x_continuous(limits = lims[[1]], expand = c(0, 0)) + 
        scale_y_continuous(limits = lims[[2]], expand = c(0, 0)) + 
        facet + ylab + gate_geom + geom_perc + 
        theme_bw() + theme(aspect.ratio = 1,
            panel.grid.minor = element_blank(), 
            axis.text = element_text(color = "black"),
            strip.background = element_rect(fill = "white"))
}
   