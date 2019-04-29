#' @rdname plotMedExprs
#' @title Plot median expressions
#' 
#' @description 
#' Plots median marker expressions across samples
#' computed on arcsinh-transformed intensities.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param k 
#'   character string. Specifies the clustering to use.
#'   If \code{facet = "antigen"}, this argument will be ignored.
#' @param facet 
#'   \code{"antigen"} or \code{"cluster_id"}. Note that 
#'   the latter requires having run \code{\link{cluster}} first.
#' @param group_by 
#'   character string. Has to appear as a column name of \code{rowData(x)}. 
#'   Specifies a factor to group samples by.
#' @param shape_by
#'   character string. Has to appear as a column name of \code{rowData(x)}. 
#'   Specifies a factor to shape by.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
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
#' 
#' # plot median expressions
#' plotMedExprs(re)
#' 
#' # run clustering
#' re <- cluster(re)
#' 
#' # plot median expressions across clusters
#' plotMedExprs(re, facet="cluster_id", k="meta8")
#' 
#' @importFrom dplyr group_by_ summarize_all
#' @importFrom reshape2 melt
# ------------------------------------------------------------------------------

setMethod(f="plotMedExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x, k="meta20", 
        facet=c("antigen", "cluster_id"), 
        group_by="condition",
        shape_by=NULL) {

        facet <- match.arg(facet)
        stopifnot(is.character(group_by), group_by %in% colnames(rowData(x)))
        
        md <- metadata(x)$experiment_info
        stopifnot(is.null(shape_by) 
            | shape_by %in% colnames(rowData(x)) 
            & is.factor(md[, shape_by]))
        if (!is.null(shape_by)) {
            n <- nlevels(md[, shape_by])
            shapes <- c(16, 17, 15, 3, 7, 8) # default shapes
            if (n > 6) {
                if (n > 18) {
                    message(paste("At most 17 shapes are currently supported", 
                        "but", n, "are required. Setting 'shape_by' to NULL."))
                    shape_by <- NULL
                } else {
                    new <- setdiff(c(seq_len(16) - 1, 18), shapes)
                    shapes <- c(shapes, sample(new, n - 6))
                }
            }
        } else {
            shapes <- NULL
        }
        
        style <- list(ylab("median expression"),
            guides(color=guide_legend(override.aes=list(alpha=1))),
            geom_point(alpha=.75, aes_string(shape=shape_by),
                position=position_jitter(width=.25, height=0)),
            scale_shape_manual(values = shapes),
            geom_boxplot(alpha=.5, width=.75, fill=NA, outlier.color=NA),
            theme_bw(), theme(
                axis.text=element_text(color="black"),
                panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="lightgrey", size=.25),
                strip.background=element_rect(fill="grey90", color=NA)))
        
        if (facet == "antigen") {
            med_exprs <- data.frame(exprs(x), sample_id=sample_ids(x)) %>% 
                group_by_(~sample_id) %>% summarize_all(funs(median))
            med_exprs <- melt(med_exprs, id.vars=c("sample_id"),
                variable.name="antigen", value.name="med_expr")
        } else {
            .check_validity_of_k(x, k)
            cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
            # compute medians across samples & clusters
            med_exprs <- data.frame(exprs(x)[, state_markers(x)], 
                sample_id=sample_ids(x), cluster_id=cluster_ids) %>% 
                group_by_(~sample_id, ~cluster_id) %>% 
                summarize_all(funs(median))
            med_exprs <- melt(med_exprs, id.vars=c("sample_id", "cluster_id"),
                variable.name="antigen", value.name="med_expr")
        }
        # add metadata information
        m <- match(med_exprs$sample_id, md$sample_id)
        med_exprs <- data.frame(med_exprs, md[m, ])
        
        if (facet == "antigen") {
            ggplot(med_exprs, 
                aes_string(x=group_by, y="med_expr", col=group_by)) +
                stat_summary(fun.y="mean", geom="point", size=2.5, shape=21) +
                facet_wrap(facets="antigen", scales="free_y") + style + theme(
                axis.text.x=element_blank(), axis.ticks.x=element_blank())
        } else {
            ggplot(med_exprs, 
                aes_string(x="antigen", y="med_expr", col=group_by)) +
                facet_wrap(facets="cluster_id", scales="free_y") + style + 
                theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
        }
    }
)
        