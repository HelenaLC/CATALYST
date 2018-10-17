#' @rdname plotClusterExprs
#' @title Plot expression distributions by cluster
#' 
#' @description 
#' Plots smoothed densities of arcsinh-transformed 
#' marker intensities by cluster.
#'
#' @param x 
#'   a \code{\link{daFrame}}.
#' @param k 
#'   character string. Specifies the clustering to use.
#' @param markers
#'   character string specifying which markers to include. Defaults to NULL 
#'   (= all markers). Alternatively, if the \code{colData(x)$marker_class} 
#'   column is specified, can be one of "type", "state", or "none".
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
#' re <- cluster(re)
#' plotClusterExprs(re)
#' 
#' @import ggplot2
#' @importFrom dplyr %>% group_by_
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom reshape2 melt
#' @importFrom S4Vectors metadata
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment colData rowData
# ------------------------------------------------------------------------------

setMethod(f="plotClusterExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x, k="meta20", markers=NULL) {
        
        # check that cluster() has been run
        md_cols <- c("SOM_codes", "cluster_codes", "delta_area")
        stopifnot(md_cols %in% names(metadata(x)))
        stopifnot("cluster_id" %in% colnames(rowData(x)))
        
        # check validity of argument 'k'
        k <- check_validity_of_k(x, k)
        
        # check validity of argument 'markers'
        if (is.null(markers)) {
            markers <- colnames(x)
        } else if (length(markers) == 1 &&
                markers %in% levels(colData(x)$marker_class)) {
            idx <- colData(x)$marker_class == markers
            if (!any(idx))
                stop(sprintf("No markers matched marker class '%s'.", markers))
            markers <- colnames(x)[idx]
        } else {
            # replace problematic characters
            markers <- gsub("-", "_", markers)
            markers <- gsub(":", ".", markers)
            stopifnot(all(markers %in% colnames(exprs(x))))
        }
        
        # calculate median markers expressions
        cluster_ids <- cluster_codes(x)[cluster_ids(x), k]
        dat <- data.frame(exprs(x)[, markers], rowData(x))
        dat$cluster_id <- cluster_ids
        meds <- dat %>% group_by_("cluster_id") %>% 
            summarise_at(markers, funs(median))
        
        # get cluster frequencies
        freqs <- table(cluster_ids) / nrow(x)
        freqs <- round(freqs * 100, 2)
        
        # constrcut data.frame & reference for plotting
        dat <- melt(dat, id.vars=names(rowData(x)), id.var="cluster_id",
            variable.name="antigen", value.name="expression")
        dat$ref <- "no"
        ref <- dat
        ref$cluster_id <- "ref"
        ref$ref <- "yes"
        df <- rbind(dat, ref)
        
        # hierarchical clustering
        d <- dist(meds, method="euclidean")
        h <- hclust(d, method="average")
        o <- h$order
        
        # reorder clusters
        df$cluster_id <- factor(df$cluster_id, 
            levels=rev(c("ref", levels(cluster_ids)[o])),
            labels=rev(c("ref", paste0(names(freqs), " (", freqs, "%)")[o])))
        
        ggplot(df, aes_string(x="expression", y="cluster_id", 
            col="ref", fill="ref")) + geom_density_ridges(alpha=.2) +
            facet_wrap(~antigen, scales="free_x", nrow=2) + 
            theme_ridges() + theme(legend.position="none",
                strip.background=element_blank(),
                strip.text=element_text(face="bold"))
    }
)
