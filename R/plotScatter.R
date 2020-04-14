#' @rdname plotScatter
#' @title Scatter plot
#' 
#' @description Bivariate scatter plots including visualization of
#' (group-specific) gates, their boundaries and percentage of selected cells.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param chs character string pecifying which channels to plot.
#'   Valid values are antigens: \code{rownames(x)}, 
#'   channel names: \code{rowData(x)$channel_name} or 
#'   other channels stored in \code{names([int_]colData(x))},
#'   and should correspond to numeric variables.
#' @param color_by character string specifying a non-numeric
#'   cell metadata column to color by; valid values are 
#'   \code{names(colData(x))}; or NULL to color by density.
#' @param facet_by character string specifying a non-numeric
#'   cell metadata column to facet by; valid values are 
#'   \code{names(colData(x))}. When \code{length(chs) == 1}, 
#'   2 facetting variables may be provided, otherwise 1 only.
#' @param bins numeric of length 1 giving the number of bins 
#'   for \code{\link[ggplot2]{geom_hex}} when coloring by density.
#' @param assay character string specifying which assay data to use.
#'   Should be one of \code{assayNames(x)}.
#' @param label character string specifying axis labels should include
#'   antigen, channel names, or a concatenation of both.
#' @param zeros logical specifying whether to include 0 values.
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(raw_data)
#' sce <- prepData(raw_data)
#' 
#' dna_chs <- c("DNA1", "DNA2")
#' plotScatter(sce, dna_chs, label = "both")
#' 
#' plotScatter(sce, 
#'   chs = sample(rownames(sce), 4), 
#'   color_by = "sample_id")
#'   
#' sce <- prepData(sample_ff)  
#' ids <- sample(rownames(sample_key), 3)
#' sce <- assignPrelim(sce, sample_key[ids, ])
#' sce <- sce[, sce$bc_id %in% ids]
#' 
#' chs <- sample(rownames(sce), 5)
#' plotScatter(sce, chs, color_by = "bc_id")
#' plotScatter(sce, chs, color_by = "delta")
#' 
#' @import ggplot2
#' @importFrom Matrix rowSums
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom stats reformulate
#' @importFrom SingleCellExperiment int_colData
#' @importFrom SummarizedExperiment assay assayNames colData
#' @export

plotScatter <- function(x, chs, color_by = NULL, facet_by = NULL,
    bins = 100, assay = "exprs", 
    label = c("antigen", "channel", "both"),
    zeros = FALSE) {
    # check validity of input arguments
    stopifnot(is.null(color_by) || is.character(color_by) && length(color_by == 1))
    #args <- as.list(environment())
    #.check_args_plotScatter(args)
   
    # subset features to speed up matrix transpose & 
    x <- x[match(chs, rownames(x), nomatch = 0), ]
    y <- assay(x, assay)[, , drop = FALSE]
    
    # rename features for visualization to
    # include both channel name & description
    rd <- rowData(x)
    m <- rd$marker_name
    c <- rd$channel_name
    nms <- switch(match.arg(label), antigen = m, channel = c,
        both = ifelse(c == m, c, paste(c, m, sep = "-")))
    rename <- chs %in% rownames(x)
    chs[rename] <- rownames(y) <- nms
    
    # construct data.frame of specified assay data & all cell metadata
    cd <- cbind(colData(x), int_colData(x))
    df <- data.frame(
        t(as.matrix(y)), cd,
        check.names = FALSE, 
        stringsAsFactors = FALSE)
    cd_vars <- intersect(names(cd), names(df))
    
    # initialize facetting & (optionally) melt data.frame 
    if (length(chs) > 2) {
        df <- melt(df, id.vars = c(chs[1], cd_vars))
        facet <- "variable"
        ylab <- ylab(NULL)
        chs[2] <- "value"
    } else {
        facet <- NULL
        ylab <- NULL
    }
    if (is.null(color_by)) {
        col_var <- guides <- NULL
        fill_var <- "..ncount.."
        scales <- scale_fill_gradientn(trans = "sqrt",
            colors = c("navy", rev(brewer.pal(11, "Spectral"))))
        geom <- geom_hex(bins = bins, na.rm = TRUE, show.legend = FALSE)
    } else {
        fill_var <- NULL
        col_var <- sprintf("`%s`", color_by)
        geom <- geom_point(alpha = 0.2, size = 0.8, na.rm = TRUE)
        if (is.numeric(df[[color_by]])) {
            guides <- NULL
            scales <- scale_color_gradientn(
                colors = c("navy", rev(brewer.pal(11, "Spectral"))))
        } else {
            scales <- NULL
            guides <- guides(col = guide_legend(
                override.aes = list(alpha = 1, size = 3)))
        }
    }
    facet <- c(facet, facet_by)
    if (!is.null(facet)) {
        if (length(facet) == 1) {
            facet <- facet_wrap(facet)    
        } else {
            facet <- facet_grid(
                cols = vars(!!sym(facet[1])), 
                rows = vars(!!sym(facet[2])))
        }
    }
    xy <- sprintf("`%s`", chs)
    if (!zeros) df <- df[rowSums(df[, chs[c(1, 2)]] == 0) == 0, ]
    ggplot(df, aes_string(xy[1], xy[2], col = col_var, fill = fill_var)) + 
        geom + scales + guides + facet + ylab + 
        theme_bw() + theme(aspect.ratio = 1,
            panel.grid = element_blank(), 
            axis.text = element_text(color = "black"),
            strip.background = element_rect(fill = "white"),
            legend.key.height = unit(0.8, "lines"))
}
