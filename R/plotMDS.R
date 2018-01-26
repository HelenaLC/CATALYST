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
#' @param color_by a character string that appears as a column name in the
#' metadata-table of the input \code{daFrame}. Specifies the color coding.
#' In the case of multiple conditions, \code{"condition"} will color points
#' according to each sample's combined conditions.
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotMDS(re)
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import ggplot2 ggrepel SummarizedExperiment
#' @importFrom dplyr group_by summarize_all
#' @importFrom ggrepel geom_label_repel
#' @importFrom limma plotMDS
# ==============================================================================

setMethod(f="plotMDS", 
    signature=signature(x="daFrame"), 
    definition=function(x, color_by="condition") {
        
        valid_cols <- colnames(rowData(x))
        if (!color_by %in% c("condition", valid_cols)) {
            txt <- dQuote(c(setdiff(valid_cols, "condition"), "condition"))
            stop("Argument 'color_by = ", dQuote(color_by), "' invalid.\n",
                "       Should be one of: ", paste(txt, collapse=", "))
        }
        tibble <- data.frame(sample_id=sample_ids(x), exprs(x)) %>% 
            group_by(sample_id) %>% summarize_all(funs(median))
        med_es <- t(tibble[, -1])
        colnames(med_es) <- tibble$sample_id
        mds <- plotMDS(med_es, plot=FALSE)
        md <- metadata(x)[[1]]
        conds <- grep("condition", colnames(md), value=TRUE)
        conds_combined <- apply(md[, conds], 1, paste, collapse="/")
        df <- data.frame(MDS1=mds$x, MDS2=mds$y, md, condition=conds_combined)
      
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
