# ==============================================================================
# Estimate separation cutoffs
# ------------------------------------------------------------------------------

#' @rdname estCutoffs
#' @title Estimation of distance separation cutoffs
#' 
#' @description 
#' For each barcode, estimates a cutoff parameter for the 
#' distance between positive and negative barcode populations.
#'
#' @param x       
#' a \code{\link{dbFrame}}.
#' 
#' @details 
#' For the estimation of sample-specific cutoff parameters, we fit a 
#' three-parameter log-logistic function to the yields function. 
#' As an adequate cutoff estimate, we target a point which approximately marks 
#' the end of the plateau regime and on-set of yield decline. 
#' By default, we compute this as the first minimum of the fifth derivative 
#' (\code{deriv=5}). However, depending on the overall doublet-singlet 
#' separation, and how bimodel the count distribution is, another derivate 
#' may provide a better estimate. 
#' As a general guideline, higher values of \code{deriv} will shift the 
#' computed minimum towards the left and yield more liberal (low) cutoffs, 
#' while a low \code{deriv} will shift it towards the inflection point of 
#' the yields functions, arriving at more stringent (high) estimates.
#'
#' @return
#' Will update the \code{sep_cutoffs} slot of the input \code{\link{dbFrame}} 
#' and return the latter.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' # assign preliminary IDs
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' # estimate separation cutoffs
#' re <- estCutoffs(x = re)
#' # view exemplary estimate
#' plotYields(re, "A1")
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom stats lm coef D
#' @importFrom drc drm LL.3
# ------------------------------------------------------------------------------

setMethod(f="estCutoffs", 
    signature=signature(x="dbFrame"), 
    definition=function(x) {
        
        sep_cutoffs <- seq(0, 1, .01)
        n_bcs <- nrow(bc_key(x))
        ests <- numeric(n_bcs)
        
        # three-parameter log-logistic function & 1st derivative
        llf <- quote(d/(1+exp(b*(log(sep_cutoffs)-log(e)))))
        deriv <-  stats::D(llf, "sep_cutoffs")
        
        for (i in seq_len(n_bcs)) {
            df <- data.frame(x=sep_cutoffs, y=as.vector(yields(x)[i, ]))
            fit <- drc::drm(y~x, data=df, fct=drc::LL.3())
            b <- fit$coefficients[1]
            d <- fit$coefficients[2]
            e <- fit$coefficients[3]
            linear_fit <- lm(yields(x)[i, ] ~ sep_cutoffs + 1)
            intercept <- coef(linear_fit)[1]
            slope <- coef(linear_fit)[2]
            rss_linear <- sum((yields(x)[i,] - predict(linear_fit)) ^ 2)
            rss_llf <- sum((yields(x)[i,] - eval(llf)) ^ 2) 
            w <- rss_llf / (rss_llf + rss_linear)
            est_llf <- sep_cutoffs[abs(c(0, eval(deriv)[-1])) / 
                    eval(llf) > 1e-2][1]
        }
        names(ests) <- rownames(bc_key(x))
        sep_cutoffs(x) <- ests
        x
    })