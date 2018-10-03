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
#'   numeric or character string. Specifies the clustering to use.
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
#' plotClusterExprs(re)
#' 
#' @import ggplot2
#' @import ggridges
#' @importFrom reshape2 melt
# ------------------------------------------------------------------------------

setMethod(f="plotClusterExprs", 
    signature=signature(x="daFrame"), 
    definition=function(x, k=20, markers=NULL) {
        
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
            stopifnot(all(markers %in% colnames(exprs(x))))
        }
        
        dat <- data.frame(exprs(x)[, markers], rowData(x))
        cluster_ids <- cluster_codes(x)[cluster_ids(x), k]
        dat$cluster_id <- cluster_ids
        dat <- melt(dat, id.vars=names(rowData(x)), id.var="cluster_id",
            variable.name="antigen", value.name="expression")
        dat$ref <- "no"
        ref <- dat
        ref$cluster_id <- "ref"
        ref$ref <- "yes"
        df <- rbind(dat, ref)
        
        # reorder clusters
        cluster_ids <- as.character(levels(cluster_ids))
        ids <- suppressWarnings(as.numeric(cluster_ids))
        ids <- c("ref", sort(cluster_ids[is.na(ids)]), sort(ids[!is.na(ids)]))
        df$cluster_id <- factor(df$cluster_id, levels=rev(ids))
        
        ggplot(df, aes_string(x="expression", y="cluster_id", 
            col="ref", fill="ref")) + geom_density_ridges(alpha=.2) +
            facet_wrap(~antigen, scales="free_x", nrow=2) + 
            theme_ridges() + theme(legend.position="none")
    }
)
