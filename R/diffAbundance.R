# ==============================================================================
# Differential cell population abundance
# ------------------------------------------------------------------------------
#' @rdname diffAbundance
#' @title Differential population abundances
#' 
#' @description 
#' Compares the proportions of cell types across experimental conditions 
#' and aims to highlight populations that are present at different ratios.
#'
#' @param x a \code{\link{daFrame}}.
#' @param k specifies which clustering to perform the analysis on.
#' @param K contrasts.
#' 
#' @return Writes a table with frequencies across samples,
#' non-adjusted and adjuster p-values for each cluster into the
#' \code{metadata} slot of the input \code{daFrame}.
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
#' data(PBMC_fs, PBMC_panel, PBMC_md)
#' re <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)
#' # specify contrasts
#' K <- matrix(c(0, 1), nrow=1, byrow=TRUE, dimnames=list("BCRXLvsRef"))
#' re <- diffAbundance(re, 12, K)
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import SummarizedExperiment
#' @importFrom lme4 glmer
#' @importFrom multcomp glht
#' @importFrom stats p.adjust
#' @export
# ==============================================================================

setMethod(f="diffAbundance", 
    signature=signature(x="daFrame"),
    definition=function(x, k, K) {

        md <- metadata(x)[[1]]
        cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
        counts <- as.data.frame.matrix(table(cluster_ids, sample_ids(x)))
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
        # get p-values
        nm <- rownames(K)
        p_vals <- do.call(rbind, fit_binomial) 
        colnames(p_vals) <- paste0("pval_", nm) 
        rownames(p_vals) <- rownames(counts)
        # get adjusted p-values
        adj_p <- apply(p_vals, 2, p.adjust, method="BH") 
        colnames(adj_p) <- paste0("adjp_", nm)
        
        re <- data.frame(cluster_id=rownames(freqs), 
            freqs, p_vals, adj_p, row.names=NULL)
        metadata(x)$diff_abundance[[nm]] <- re
        return(x)
    }
)
