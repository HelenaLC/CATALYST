#' @rdname plotSNE
#' @title plot t-SNE
#' 
#' @description t-SNE plot colored by marker expression or clustering. 
#'
#' @param x a \code{\link{daFrame}}.
#' @param color_by numeric value or character string specifying a clustering
#' OR a character string specifying an antibody whose expression to color by.
#' @param facet a character string specifying a factor to subset the data by.
#' One of \code{names(rowData(x))}. Defaults to NULL.
#' 
#' @return a \code{ggplot} object.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' # construct daFrame
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # run clustering
#' lineage <- c("CD3", "CD45", "CD4", "CD20", "CD33", 
#'     "CD123", "CD14", "IgM", "HLA_DR", "CD7")
#' re <- cluster(re, cols_to_use=lineage)
#' 
#' # run t-SNE
#' re <- tSNE(re, n=100)
#' 
#' # color by clustering
#' plotSNE(re, color_by=12)
#' 
#' # color by marker expression
#' plotSNE(re, color_by="pS6", facet="condition")
#' 
#' @import ggplot2
# ------------------------------------------------------------------------------

setMethod(f="plotSNE",
    signature=signature(x="daFrame"),
    definition=function(x, color_by=20, facet=NULL) {
        
        # check validity of 'color_by' argument
        color_by <- as.character(color_by)
        if (!color_by %in% colnames(exprs(x)))
            check_validity_of_k(x, color_by)
        # check validity of 'facet' argument
        if (!is.null(facet) && !facet %in% names(rowData(x)))
            stop("'facet' should be either NULL, or one of\n",
                paste(dQuote(colnames(rowData(x))), collapse=", "))
        
        tsne_inds <- metadata(x)$tsne_inds
        tsne <- metadata(x)$tsne
        
        df <- data.frame(
            exprs(x)[tsne_inds, ], 
            rowData(x)[tsne_inds, ],
            tSNE1=tsne$Y[, 1], tSNE2=tsne$Y[, 2])

        p <- ggplot(df, aes_string(x="tSNE1", y="tSNE2")) +
            theme_void() + theme(aspect.ratio=1,
                panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="grey", size=.25),
                axis.text=element_text(color="black"),
                axis.title=element_text(color="black", face="bold"))
        
        if (color_by %in% colnames(exprs(x))) {
            pal <- rev(brewer.pal(11, "Spectral"))
            p <- p + geom_point(size=.75, aes_string(color=color_by)) +
                scale_color_gradientn(color_by, colors=pal)
        } else {
            cluster_ids <- cluster_codes(x)[, color_by][cluster_ids(x)[tsne_inds]]
            df$cluster_id <- factor(cluster_ids)
            cols <- cluster_cols[seq_along(unique(df$cluster_id))]
            names(cols) <- levels(cols)
            if (length(cols) > 10) n_col <- 2 else n_col <- 1
            p <- p + geom_point(data=df, size=.75, 
                aes_string(color="cluster_id")) +
                guides(color=guide_legend(override.aes=list(size=3),
                    ncol=n_col)) + scale_color_manual(values=cols)
        }
        if (is.null(facet)) {
            p
        } else {
            p + facet_wrap(facet)
        }
    }
)
