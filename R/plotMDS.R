# ==============================================================================
# MDS plot on median marker expression values
# ------------------------------------------------------------------------------

#' @rdname plotMDS
#' @title MDS plot
#' 
#' @description 
#' Multi-dimensional scaling (MDS) plot on median marker expressions.
#'
#' @param x a \code{\link{daFrame}}.
#' 
#' @return 
#' Plots a multi-dimensional scaling (MDS) plot 
#' on median marker expression values.
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
#' @importFrom flowCore flowSet fsApply
#' @importFrom limma plotMDS
# ==============================================================================

setMethod(f="plotMDS", 
    signature=signature(x="daFrame"), 
    definition=function(x) {
        # get parameters descriptions, replace dash by underscore
        d <- parameters(data(x)[[1]])$desc
        d <- gsub("-", "_", d)
        # get expressions, use antigens as column names
        es <- fsApply(data(x), exprs)
        colnames(es) <- d
        # extract lineage & functional markers
        inds <- panel(x)$Antigen[panel(x)$Lineage | panel(x)$Functional]
        # arcsinh transformation & column subsetting
        es <- asinh(es[, inds]/5)
        # compute median expressions
        n_events <- fsApply(data(x), nrow)
        df <- data.frame(sample_id=rep(metadata(x)$sample_id, n_events), es)
        tibble <- df %>% group_by(sample_id) %>% summarize_all(funs(median))
        med_es <- t(tibble[, -1])
        colnames(med_es) <- tibble$sample_id
        mds <- limma::plotMDS(med_es, plot=FALSE)
        df <- data.frame(MDS1=mds$x, MDS2=mds$y, 
            sample_id=metadata(x)$sample_id, 
            condition=metadata(x)$condition)

        r <- diff(range(df$MDS2)) / diff(range(df$MDS1))
        ggplot(df, aes_string(x="MDS1", y="MDS2", col="condition")) + 
            geom_point(size=10, alpha=.75) + 
            geom_label_repel(aes_string(label="sample_id"), show.legend=FALSE) + 
            theme_void() + theme(aspect.ratio=r,
                panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color='lightgrey', size=.25), 
                axis.title=element_text(face='bold'),
                axis.text=element_text())
    }
)