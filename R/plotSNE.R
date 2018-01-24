# ==============================================================================
# tSNE-method for dbFrame-class
# ------------------------------------------------------------------------------
#' @rdname plotSNE
#' @title MDS plot
#' 
#' @description 
#' Multi-dimensional scaling (MDS) plot on median marker expressions.
#'
#' @param x a \code{\link{daFrame}}.
#' @param color numeric value between 2 and 20 OR a character string specifying
#' an antibody that appears in the metadata table of the input \code{daFrame}.
#' Specifies the color coding.
#' @param facette one of \code{NULL}, \code{"sample_id"} or \code{"condition"}.
#' Specifies whether the data should be split between sample IDs or conditions,
#' respectively. Defaults to NULL.
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- tSNE(re)
#' plotSNE(re, "CD4")
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import ggplot2
# ==============================================================================

setMethod(f="plotSNE",
    signature=signature(x="daFrame"),
    definition=function(x, color=20, facette=NULL) {
        # check validity of 'color' argument
        color <- as.character(color)
        if (!color %in% c(colnames(exprs(x)),
            colnames(metadata(x)$cluster_codes))) 
            stop("'color' should be either a numeric value between 2 and 20,
                a character string specifying a merging to use,\n", 
                "or a character string that corresponds to a lineage marker.")
        # check validity of 'facette' argument
        if (!is.null(facette) && !facette %in% c("sample", "condition"))
            stop("'facette' should be either NULL,",
                " \"sample\" or \"condition\".")
        
        tsne_inds <- metadata(x)$tsne_inds
        tsne <- metadata(x)$tsne
        
        df <- data.frame(
            exprs(x)[tsne_inds, ],
            tSNE1=tsne$Y[, 1], tSNE2=tsne$Y[, 2],
            sample_id=sample_ids(x)[tsne_inds],
            condition=rowData(x)$condition[tsne_inds])
        p <- ggplot(df, aes_string(x="tSNE1", y="tSNE2")) +
            theme_void() + theme(aspect.ratio=1,
                panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="grey", size=.25),
                axis.text=element_text(color="black"),
                axis.title=element_text(color="black", face="bold"))
        
        if (color %in% colnames(exprs(x))) {
            pal <- rev(brewer.pal(11, "Spectral"))
            p <- p + geom_point(size=.75, aes_string(color=color)) +
                scale_color_gradientn(color, colors=pal)
        } else {
            cluster_ids <- cluster_ids(x)[tsne_inds]
            df$cluster_id <- factor(cluster_codes(x)[, color][cluster_ids])
            cols <- cluster_cols[seq_len(nlevels(df$cluster_id))]
            names(cols) <- levels(cols)
            if (length(cols) > 10) n_col <- 2 else n_col <- 1
            p <- p + geom_point(data=df, size=.75, 
                aes_string(color="cluster_id")) +
                guides(color=guide_legend(override.aes=list(size=3),
                    ncol=n_col)) + scale_color_manual(values=cols)
        }
        if (is.null(facette)) {
            p
        } else if (facette == "sample") {
            p + facet_wrap(~sample_id)
        } else if (facette == "condition") {
            p + facet_wrap(~condition)
        }
    }
)
