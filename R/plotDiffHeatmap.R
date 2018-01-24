# ==============================================================================
# Heatmap of median marker expressions across samples
# ------------------------------------------------------------------------------
#' @rdname plotDiffHeatmap
#' @title Median marker expressions across samples
#'
#' @param x a \code{\link{daFrame}}.
#' @param type a character string. Specifies the type of analysis:
#' differential \code{"abundance"} or \code{"expr"}.
#' @param contrast a character string that specifies which contrast to plot.
#' @param th a numeric giving the threshold for Z-normalization.
#' @param FDR_cutoff a numeric giving the FDR cutoff.
#' @param main a character string to use as the plot title.
#' @param palette a character string specifying the colors to interpolate.
#' @param out_path a character string. If specified, 
#' output will be generated in this location. Defaults to NULL.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @seealso \code{\link{diffAbundance}}, \code{\link{diffExpr}}
#' 
#' @examples
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' 
#' # specify contrasts
#' K <- matrix(c(0, 1), nrow=1, byrow=TRUE, dimnames=list("BCRXLvsRef"))
#' 
#' # differential abundances
#' re <- diffAbundance(re, 12, K)
#' plotDiffHeatmap(re, type="abundance", contrast="BCRXLvsRef")
#' 
#' # differential expression
#' re <- diffExpr(re, 12, K)
#' plotDiffHeatmap(re, type="expr", contrast="BCRXLvsRef", FDR_cutoff=.5)
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import pheatmap SummarizedExperiment
#' @importFrom dplyr funs group_by summarize_all
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats hclust sd
#' @export
# ==============================================================================

setMethod(f="plotDiffHeatmap", 
    signature=signature(x="daFrame"), 
    definition=function(x, type=c("abundance", "expr"), 
        contrast, clustering, th=2.5, FDR_cutoff=.05, main=NULL, 
        palette=brewer.pal(n=11, name="RdYlGn"), out_path=NULL) {
        
        # get differential analysis results
        type <- match.arg(type)
        da <- metadata(x)[[paste0("diff_", type)]]
        if (is.null(da)) {
            s <- strsplit(type, " ")[[1]]
            s <- paste0(toupper(substring(s,1,1)), substring(s,2))
            stop("No differential analysis of type '", type, "' available.", 
                "\n       Please run 'diff", s, "()' first.")
        }
        if (is.null(da[[contrast]]))
            stop("Differential ", type, " analysis '", contrast, 
                "' doesn't exist.", "\n       Available analyses are: ", 
                paste(names(da), collapse=", "))
        da <- da[[contrast]]
        
        # get significant clusters / markers & sort
        adj_p <- setNames(da[[paste0("adjp_", contrast)]], rownames(da))
        s <- names(which(sort(adj_p) < FDR_cutoff))
        if (length(s) == 0) {
            txt <- switch(type, 
                abundance = "abundant populations", 
                expr = "expressed markers")
            stop("No significantly differentially ", txt,
                " detected for 'FDR_cutoff = ", FDR_cutoff, "'.")
        }

        data <- switch(type,
            abundance = {
                # get cluster IDs
                cluster_ids <- cluster_codes(x)[, clustering][cluster_ids(x)]
                # compute cluster frequencies across samples
                counts <- as.data.frame.matrix(
                    table(cluster_ids, sample_ids(x)))
                freqs <- t(t(counts) / colSums(counts))
                # arcsine-square-root transformation
                freqs_t <- as.data.frame.matrix(asin(sqrt(freqs)))
                # Z-score normalization
                Z_norm(freqs_t[s, ], th)
            },
            expr = {
                # compute medians across samples & clusters
                md <- metadata(x)[[1]]
                cluster_ids <- cluster_codes(x)[, clustering][cluster_ids(x)]
                med_exprs <- data.frame(exprs(x), 
                    sample_id=sample_ids(x), cluster_id=cluster_ids) %>% 
                    group_by(sample_id, cluster_id) %>% 
                    summarize_all(funs(median))
                med_exprs <- melt(med_exprs, id.vars=c("sample_id", "cluster_id"),
                    variable.name="antigen", value.name="med_expr")
                med_exprs <- dcast(med_exprs, value.var = "med_expr", 
                    formula=cluster_id+antigen~sample_id)
                rownames(med_exprs) <- paste0(
                    med_exprs$cluster_id, "_", med_exprs$antigen)

                # group by cluster
                adj_p <- adj_p[s]
                o <- order(med_exprs[s, "cluster_id"], adj_p)
                s <- s[o]
                # Z-score normalization
                med_exprs_s <- med_exprs[s, levels(metadata(x)[[1]]$sample_id)]
                Z_norm(med_exprs_s, th)
            })
        
        # row labels
        row_labs <- paste0(s, " (", sprintf("%.02e", adj_p), ")")

        if (is.null(main)) 
            main <- switch(type, 
                abundance = "Normalized population proportions",
                expr = "Normalized functional marker expressions")
        
        pheatmap(data, main=main, border_color="white",
            gaps_col=as.numeric(table(metadata(x)[[1]]$condition)),
            color=colorRampPalette(palette)(100),
            breaks=seq(-th, th, length.out=101), legend_breaks=(-th):th,
            labels_row=row_labs, cluster_cols=FALSE, cluster_rows=FALSE)
    }
)
