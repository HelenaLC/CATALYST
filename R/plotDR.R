#' @rdname plotDR
#' @title Plot dim. reduction from a \code{daFrame}
#' 
#' @description Plot cell-level reduced dimensions
#'   stored within a \code{daFrame} object.
#' 
#' @param x a \code{\link{daFrame}}.
#' @param dr character string specifying the dimensionaly reduction method.
#' @param color_by character string specifying a clustering,
#'   marker, or \code{rowData} column name to color by.
#' @param facet a character string specifying 
#'   a \code{rowData} column to facet by.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' daf <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' daf <- runDR(daf, "PCA")
#' plotDR(daf, "PCA", color_by = "condition")
#' 
#' daf <- cluster(daf)
#' plotDR(daf, "PCA", color_by = "meta5")
#' 
#' @return a \code{ggplot} object.
#'   
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch}
#'  
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom SingleCellExperiment int_elementMetadata reducedDimNames
#' @importFrom SummarizedExperiment rowData

setMethod("plotDR",
    signature=signature(x = "daFrame"),
    definition=function(x, 
        dr=c("TSNE", "PCA", "MDS", "UMAP", "DiffusionMap"), 
        color_by = "meta20", facet = NULL) {
        
        dr <- match.arg(dr)
        stopifnot(dr %in% reducedDimNames(x))
        
        # check validity of 'color_by' argument
        if(!color_by %in% c(colnames(x), 
            names(cluster_codes(x)), names(rowData(x))))
            stop("Argument 'color_by' invalid. Should correspond\n", 
                "to either a marker, clustering, or rowData column.")
        if (color_by %in% names(cluster_codes(x)))
            .check_validity_of_k(x, color_by)

        # check validity of 'facet' argument
        if (!is.null(facet) && !facet %in% names(rowData(x)))
            stop("'facet' should be either NULL, or one of\n",
                paste(dQuote(names(rowData(x))), collapse=", "))
        
        # get reduced dimensions & cell indices used for computation
        idx <- int_elementMetadata(x)
        idx <- idx[[sprintf("idx.%s", dr)]]
        x <- x[unlist(idx), ]
        dr <- reducedDim(x, dr)
        df <- data.frame(dr, exprs(x), rowData(x), check.names = FALSE)
        
        p <- ggplot(df, aes_string(x=colnames(dr)[1], y=colnames(dr)[2])) +
            theme_void() + theme(aspect.ratio=1,
                panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="grey", size=.2),
                axis.text=element_text(color="black"),
                axis.title=element_text(color="black", face="bold"))
        
        if (color_by %in% colnames(x)) {
            p <- p +
                geom_point(size=.8, aes_string(color=color_by)) +
                scale_color_viridis_c()
        } else if (color_by %in% names(cluster_codes(x))) {
            # get cluster IDs
            df$cluster_id  <- .get_cluster_ids(x, color_by)
            # expand palette if more than 30 clusters
            nk <- nlevels(df$cluster_id )
            if (nk > 30) {
                cols <- colorRampPalette(.cluster_cols)(nk)
            } else {
                cols <- .cluster_cols[seq_len(nk)]
            }
            names(cols) <- levels(cluster_ids)
            if (nk > 10) n_col <- 2 else n_col <- 1
            p <- p + 
                geom_point(data=df, size=.8, aes_string(color="cluster_id")) +
                guides(color=guide_legend(override.aes=list(size=3),
                    ncol=n_col)) + scale_color_manual(values=cols)
        } else (
            p <- p + 
                geom_point(data=df, size=.8, aes_string(color=color_by)) +
                guides(color=guide_legend(override.aes=list(size=3)))
        )
        if (is.null(facet)) p else p + facet_wrap(facet, ncol=4)
    }
)
