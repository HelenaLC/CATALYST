#' @rdname plotCounts
#' @title Plot cell counts
#' 
#' @description Barplot of the number of cells measured for each sample.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param color_by 
#'   character string. Must appear as a column name of \code{rowData(x)}. 
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
#' @return a \code{\link{ggplot}} object.
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' plotCounts(re)
#' 
#' @import ggplot2 SummarizedExperiment
# ------------------------------------------------------------------------------

setMethod(f="plotCounts", 
    signature=signature(x="daFrame"), 
    definition=function(x, color_by="condition") {
        
        md <- metadata(x)$experiment_info
        valid <- names(md)
        if (!color_by %in% valid)
            stop("Argument 'color_by = ", dQuote(color_by), "' invalid.\n",
                "Should be one of: ", paste(dQuote(valid), collapse=", "))
        
        max <- ceiling(max(n_cells(x))/1e3)*1e3 + 1e3
        df <- data.frame(n_cells(x), md)
        
        ggplot(df, aes_string(x="sample_id", y="n_cells",
            fill=color_by)) + geom_bar(stat="identity", width=.75) +  
            geom_label(aes_string(label="n_cells", vjust=-.1), 
                fontface="bold.italic", fill=NA, label.size=0) +
            scale_y_continuous(limits=c(0,max), expand=c(0,0)) +
            theme_minimal() + theme(
                panel.grid.minor=element_blank(), 
                panel.grid.major.x=element_blank(),
                panel.grid.major.y=element_line(color="grey", size=.3), 
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, hjust=1, vjust=1))
    }
)
