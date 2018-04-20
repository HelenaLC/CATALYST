#' @rdname plotMDS
#' @title MDS plot
#' 
#' @description 
#' Multi-dimensional scaling (MDS) plot on median marker expressions.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param color_by 
#'   character string that appears as a column name of \code{rowData(x)}.
#'   Specifies the color coding.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotMDS(re)
#' 
#' @import ggplot2 ggrepel SummarizedExperiment
#' @importFrom dplyr group_by_ summarize_all
#' @importFrom ggrepel geom_label_repel
#' @importFrom limma plotMDS
#' @importFrom magrittr %>%
# ------------------------------------------------------------------------------

setMethod(f="plotMDS", 
    signature=signature(x="daFrame"), 
    definition=function(x, color_by="condition") {
        valid <- colnames(rowData(x))
        if (!color_by %in% valid)
            stop("Argument 'color_by = ", dQuote(color_by), "' invalid.\n",
                "Should be one of: ", paste(dQuote(valid), collapse=", "))
        
        # compute medians across samples
        meds <- data.frame(sample_id=sample_ids(x), exprs(x)) %>% 
            group_by_(~sample_id) %>% summarize_all(funs(median))
        meds <- t(data.frame(meds, row.names=1))
        
        # get MDS coordinates
        mds <- plotMDS(meds, plot=FALSE)
        df <- data.frame(MDS1=mds$x, MDS2=mds$y)
        
        # add metadata information
        md <- metadata(x)$experiment_info
        m <- match(rownames(df), md$sample_id)
        df <- data.frame(df, md[m, ])
        
        ggplot(df, aes_string(x="MDS1", y="MDS2", col=color_by)) + 
            geom_label_repel(aes_string(label="sample_id"), 
                show.legend=FALSE) + geom_point(alpha=.75) + 
            guides(col=guide_legend(overide.aes=list(alpha=1))) +
            theme_void() + theme(aspect.ratio=1,
                panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color='lightgrey', size=.25), 
                axis.title=element_text(face='bold'),
                axis.text=element_text())
    }
)
