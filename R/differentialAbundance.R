# ==============================================================================
# Differential cell population abundance
# ------------------------------------------------------------------------------
#' @rdname differentialAbundance
#' @title differentialAbundance
#' 
#' @description 
#' Compares the proportions of cell types across experimental conditions 
#' and aims to highlight populations that are present at different ratios.
#'
#' @param fs a \code{\link{flowSet}} holding all samples.
#' @param md a data.frame containing the metadata.
#' @param cluster_ids a numeric vector.
#' @param K contrasts.
#' 
#' @return Returns a table with frequencies across samples,
#' non-adjusted and adjuster p-values for each cluster.
#' 
#' @details 
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
#' 
#' @importFrom lme4 glmer
#' @importFrom multcomp glht
#' @export
# ==============================================================================

setMethod(f="differentialAbundance", 
    signature=signature(fs="flowSet", md="data.frame", 
        cluster_ids="numeric", K="matrix"),
    definition=function(fs, md, cluster_ids, K) {

        sample_ids <- rep(md$sample_id, fsApply(fs, nrow))
        counts <- as.data.frame.matrix(table(mc_ids, sample_ids))
        freqs <- t(t(counts) / colSums(counts))
        n_events <- colSums(counts)
        
        # fit GLMM for ea. cluster
        n_clusters <- nrow(counts)
        formula <- counts/total ~ condition + (1|sample_id)
        fit_binomial <- lapply(seq_len(n_clusters), function(i) {
            data <- data.frame(
                counts=as.numeric(counts[i, md$sample_id]),
                total=n_events[md$sample_id], md)
            fit <- glmer(formula, weights=total, family=binomial, data=data)
            # fit contrasts one by one
            p_val <- apply(K, 1, function(k) {
                contr <- glht(fit, linfct=matrix(k, 1)) 
                summary(contr)$test$pvalues
            })
            return(p_val)
        })
        p_vals <- do.call(rbind, fit_binomial) 
        colnames(p_vals) <- paste0("pval_", contrast_names) 
        rownames(p_vals) <- rownames(counts)
        # ddjust the p-values
        adj_p <- apply(p_vals, 2, p.adjust, method="BH") 
        colnames(adj_p) <- paste0("adjp_", contrast_names)
        data.frame(cluster_id=rownames(freqs), 
            freqs, p_vals, adj_p, row.names=NULL)
    }
)