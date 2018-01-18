# ==============================================================================
# Marker ranking based on the non-redundancy score
# ------------------------------------------------------------------------------
#' @rdname plotNRS
#' @title Non-redundancy scores across markers
#' 
#' @description
#'
#' @param x a \code{\link{daFrame}}
#' 
#' @return
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
#' @importFrom matrixStats colQuantiles
#' @importFrom reshape2 melt
#' @export
# ==============================================================================

setMethod(f="plotNRS", 
    signature=signature(x="daFrame"), 
    definition=function(x) {
        
        # get parameters descriptions, replace dash by underscore
        d <- parameters(data(x)[[1]])$desc
        d <- gsub("-", "_", d)
        lineage <- panel(x)$Antigen[as.logical(panel(x)$Lineage)]
        lineage <- gsub("-", "_", lineage)
        
        # arcsinh-transformation & column subsetting
        fs <- fsApply(fs, function(ff) {
            colnames(ff) <- d
            exprs(ff) <- asinh(exprs(ff[, lineage])/5)
            return(ff)
        })
        
        # calculate NRS
        scores <- fsApply(fs[, lineage], nrs, use.exprs=TRUE)
        rownames(scores) <- metadata(x)$sample_id
        mean_scores <- colMeans(scores, na.rm=TRUE)
        
        # plot NRS in decreasing order
        o <- names(sort(mean_scores, decreasing=TRUE))
        scores <- data.frame(scores)
        scores$sample_id <- metadata(x)$sample_id
        
        df <- reshape2::melt(scores, id.var="sample_id",
            value.name="NRS", variable.name="antigen")
        df$antigen <- factor(df$antigen, levels=o)
        m <- match(df$sample_id, metadata(x)$sample_id)
        df$condition <- md$condition[m]
        
        ggplot(df, aes_string(x="antigen", y="NRS")) +
            geom_point(aes_string(color="condition"), size=2.5, alpha=.75,
                position=position_jitter(width=.25, height=0)) +
            guides(color=guide_legend(override.aes=list(alpha=1))) +
            stat_summary(fun.y="mean", geom="point", fill="white", shape=21) +
            geom_boxplot(width=.75, fill=NA, outlier.color=NA) +
            theme_classic() + theme(
                panel.grid.major=element_line(color="lightgrey", size=.25),
                axis.text=element_text(color="black"),
                axis.text.x=element_text(angle=90, vjust=.5, hjust=1))
    }
)