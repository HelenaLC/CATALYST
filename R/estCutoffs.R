#' @rdname estCutoffs
#' @title Estimation of distance separation cutoffs
#' 
#' @description 
#' For each sample, estimates a cutoff parameter for the 
#' distance between positive and negative barcode populations.
#'
#' @param x       
#'   a \code{\link{dbFrame}}.
#' 
#' @details 
#' For the estimation of cutoff parameters, we considered yields
#' upon debarcoding as a function of the applied cutoffs. 
#' Commonly, this function will be characterized by an initial weak decline, 
#' where doublets are excluded, and subsequent rapid decline in yields to zero. 
#' In between, low numbers of counts with intermediate barcode separation 
#' give rise to a plateau. As an adequate cutoff estimate, 
#' we target the point that approximately marks the end of the plateau regime 
#' and the onset of yield decline. To facilitate robust cutoff estimation, 
#' we fit a linear and a three-parameter log-logistic function 
#' to the yields function:
#' \deqn{f(x)=\frac{d}{1+e^{b(log(x)-log(e))}}}{
#' f(x) = d / (1 + exp(b * (log(x) - log(e))))}
#' The goodness of the linear fit relative to the log-logistic fit 
#' is weighed with:
#' \deqn{w=\frac{RSS_{log-logistic}}{RSS_{log-logistic}+RSS_{linear}}}{
#' w = RSS(log-logistic) / (RSS(log-logistic) + RSS(linear))}
#' and the cutoffs for both functions are defined as:
#' \deqn{c_{linear}=-\frac{\beta_0}{2\beta_1}}{
#' c(linear) = - beta0 / (2 * beta1)}
#' \deqn{c_{log-logistic}=argmin_x\{\frac{\vert f'(x)\vert}{
#' f(x)}>0.1\}}{c(log-logistic) = argmin x { | f'(x) | / f(x) > 0.1 }}
#' The final cutoff estimate is defined as the weighted mean 
#' between these estimates:
#' \deqn{c=(1-w)\cdot c_{log-logistic}+w\cdot c_{linear}}{
#' c = (1 - w) x c(log-logistic) + w x c(linear)}
#' 
#' @return
#' Will update the \code{sep_cutoffs} slot of the input \code{\link{dbFrame}} 
#' and return the latter.
#' 
#' @author Helena Lucia Crowell \email{helena.crowell@uzh.ch}
#'
#' @references 
#' Finney, D.J. (1971). Probit Analsis. 
#' \emph{Journal of Pharmaceutical Sciences} \bold{60}, 1432. 
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
            fit <- tryCatch(
                drc::drm(y~x, data=df, fct=drc::LL.3()), 
                error=function(e) e)
            if (inherits(fit, "error")) 
                next
            b <- fit$coefficients[1]
            d <- fit$coefficients[2]
            e <- fit$coefficients[3]
            linear_fit <- lm(yields(x)[i, ] ~ sep_cutoffs + 1)
            intercept <- coef(linear_fit)[1]
            slope <- coef(linear_fit)[2]
            rss_linear <- sum((yields(x)[i,] - predict(linear_fit)) ^ 2)
            rss_llf <- sum((yields(x)[i,] - eval(llf)) ^ 2) 
            w <- rss_llf / (rss_llf + rss_linear)
            est_linear <- round(-intercept/(2*slope), 2)
            est_llf <- sep_cutoffs[abs(c(0, eval(deriv)[-1])) / 
                    eval(llf) > 1e-2][1]
            ests[i] <- (1 - w) * est_llf + w * est_linear
        }
        ests <- round(ests, 2)
        names(ests) <- rownames(bc_key(x))
        sep_cutoffs(x) <- ests
        x
    })
            