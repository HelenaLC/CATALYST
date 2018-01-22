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
#' @param col numeric value between 2 and 20 OR a character string specifying
#' an antibody that appears in the metadata table of the input \code{daFrame}.
#' Specifies the color coding.
#' @param facette one of \code{NULL}, \code{"sample_id"} or \code{"condition"}.
#' Specifies whether the data should be split between sample IDs or conditions,
#' respectively. Defaults to NULL.
#' 
#' @return a \code{ggplot} object.
#'
#' @references 
#' Nowicka M, Krieg C, Weber LM et al.
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @examples
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2
# ==============================================================================

setMethod(f="plotSNE",
    signature=signature(x="daFrame"),
    definition=function(x, col=20, facette=NULL, which=100) {
        # check validity of input arguments
        invalid <- FALSE
        if (is.character(col)) {
            if (!col %in% colnames(exprs(x)))
                invalid <- TRUE
        } else if (is.numeric(col)) {
            if (as.integer(col) != col | col < 2 | col > 20) 
                invalid <- TRUE
        } 
        if (invalid)  
            stop("'col' should be either a numeric value between 2 and 20\n", 
                "or a character string that corresponds to a lineage marker.")
        if (!is.null(facette) && !facette %in% c("sample", "condition"))
            stop("'facette' should be either NULL,",
                " \"sample\" or \"condition\".")
        
        tsne_inds <- metadata(x)$tsne_inds
        tsne <- metadata(x)$tsne
        
        df <- data.frame(
            exprs(x)[tsne_inds, ],
            tSNE1=tsne$Y[, 1], tSNE2=tsne$Y[, 2],
            sample_id=sample_ids(x)[tsne_inds],
            condition=conditions(x)[tsne_inds])
        p <- ggplot(df, aes_string(x="tSNE1", y="tSNE2")) +
            theme_void() + theme(aspect.ratio=1,
                panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="grey", size=.25),
                axis.text=element_text(color="black"),
                axis.title=element_text(color="black", face="bold"))
        
        if (is.character(col)) {
            pal <- rev(brewer.pal(11, "Spectral"))
            p <- p + geom_point(size=.75, aes_string(col=col)) +
                scale_color_gradientn(col, colors=pal)
        } else if (is.numeric(col)) {
            col_nm <- paste0("k", col)
            df$cluster_id <- factor(cluster_ids(x)[tsne_inds, col_nm])
            if (col > 10) n_col <- 2 else n_col <- 1
            p <- p + geom_point(data=df, size=.75, 
                aes_string(col="cluster_id")) +
                guides(color=guide_legend(override.aes=list(size=3),
                    ncol=n_col)) + scale_color_manual(values=cluster_cols)
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