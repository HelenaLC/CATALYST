# ==============================================================================
# Differential analysis of marker expressions
# ------------------------------------------------------------------------------
#' @rdname diffExpr
#' @title Differential marker expressions
#' 
#' @description 
#' Compares the proportions of cell types across experimental conditions 
#' and aims to highlight populations that are present at different ratios.
#'
#' @param x a \code{\link{daFrame}}.
#' @param k specifies which clustering to perform the analysis on.
#' @param K contrasts.
#' @param model a character string specifying which model to use;
#' \code{"lm"} (linear model) or \code{"lmm"} (lineager mixed model)
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
#' diffExpr(re, 20, K, "lm")
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import SummarizedExperiment
#' @importFrom lme4 glmer lmer
#' @importFrom dplyr group_by summarize_all
#' @importFrom multcomp glht
#' @importFrom reshape2 dcast melt
#' @importFrom stats lm p.adjust
# ==============================================================================

setMethod(f="diffExpr", 
    signature=signature(x="daFrame"),
    definition=function(x, k, K, model=c("lm", "lmm")) {

        # compute medians across samples & clusters
        md <- metadata(x)[[1]]
        cluster_ids <- factor(cluster_codes(x)[, k][cluster_ids(x)])
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
        
        # fit LM or LMM for ea. marker
        fml_lm  <- y ~ condition
        fml_lmm <- y ~ condition + (1 | patient_id) 
        fit_gaussian <- lapply(seq_len(nrow(med_exprs)), function(i) {
            data <- data.frame(y=as.numeric(med_exprs[i, md$sample_id]), md)
            fit <- switch(match.arg(model, c("lm", "lmm")),
                lm = lm(fml_lm, data),
                lmm = lmer(fml_lmm, data))
            # fit contrasts
            apply(K, 1, function(k) {
                contrast <- glht(fit, linfct=matrix(k, 1))
                p_val <- summary(contrast)$test$pvalues
                return(p_val)
            })
        })
        # get p-values
        nm <- rownames(K)
        p_vals <- do.call(rbind, fit_gaussian)
        colnames(p_vals) <- paste0("pval_", nm)
        rownames(p_vals) <- rownames(med_exprs)
        # get adjusted p-values
        adj_p <- apply(p_vals, 2, p.adjust, method="BH")
        colnames(adj_p) <- paste0("adjp_", nm)
        metadata(x)$diff_expr[[nm]] <- data.frame(p_vals, adj_p)
        return(x)
    }
)
