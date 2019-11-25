#' @rdname estCutoffs
#' @title Estimation of distance separation cutoffs
#' @description For each sample, estimates a cutoff parameter for 
#' the distance between positive and negative barcode populations.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
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
#' @return the input \code{SingleCellExperiment} \code{x} is returned
#' with an additional \code{int_metadata} slot \code{sep_cutoffs}
#' stored in \code{altExp(x, altExp)}. If \code{altExp} is NULL, 
#' \code{sep_cutoffs} are stored in the \code{int_metadata} of \code{x}.
#' 
#' @author Helena L. Crowell
#'
#' @references Finney, D.J. (1971). Probit Analsis. 
#' \emph{Journal of Pharmaceutical Sciences} \bold{60}, 1432. 
#' 
#' @examples
#' library(SingleCellExperiment)
#' 
#' # construct SCE
#' data(sample_ff, sample_key)
#' sce <- fcs2sce(sample_ff)
#' sce <- transform(sce)
#'     
#' # assign preliminary barcode IDs
#' # & estimate separation cutoffs
#' sce <- assignPrelim(x = sce, bc_key = sample_key)
#' sce <- estCutoffs(sce)
#' 
#' # access separation cutoff estimates
#' metadata(sce)$sep_cutoffs
#' 
#' # view yield plots including current cutoff
#' plotYields(sce, which = "A1")
#' 
#' @importFrom drc drm LL.3
#' @importFrom Matrix colMeans
#' @importFrom methods is
#' @importFrom stats coef D predict
#' @importFrom S4Vectors metadata metadata<-
#' @export

estCutoffs <- function(x) {
    stopifnot(is(x, "SingleCellExperiment"),
        !is.null(metadata(x)$bc_key),
        !is.null(x$bc_id), !is.null(x$delta))
    
    ids <- rownames(bc_key <- metadata(x)$bc_key)
    n_seps <- length(names(seps) <- seps <- seq(0, 1, 0.01))
    
    # split cell by barcode ID
    cs <- split(seq_len(ncol(x)), x$bc_id)
    
    # compute yields upon applying separation cutoffs
    ys <- vapply(seps, function(u) x$delta >= u, numeric(ncol(x)))
    ys <- vapply(ids, function(id) colMeans(ys[cs[[id]], ]), numeric(n_seps))
    
    # three-parameter log-logistic function & 1st derivative
    dll <- D(ll <- quote(d/(1+exp(b*(log(seps)-log(e))))), "seps")
    
    ests <- vapply(ids, function(id) {
        df <- data.frame(x = seps, y = ys[, id])
        fit <- tryCatch(
            drm(y~x, data = df, fct = LL.3()),
            error = function(e) e)
        if (inherits(fit, "error")) 
            return(NA)
        b <- fit$coefficients[1]
        d <- fit$coefficients[2]
        e <- fit$coefficients[3]
        lm_fit <- lm(ys[, id] ~ seps + 1)
        rss_lm <- sum((ys[, id] - predict(lm_fit)) ^ 2)
        rss_ll <- sum((ys[, id] - eval(ll)) ^ 2)
        est_lm <- - coef(lm_fit)[1] / (2 * coef(lm_fit)[2])
        est_ll <- seps[abs(c(0, eval(dll)[-1])) / eval(ll) > 1e-2][1]
        w <- rss_ll / (rss_ll + rss_lm)
        w * est_lm + (1 - w) * est_ll
    }, numeric(1))
    
    # store estimates in metadata & return SCE
    metadata(x)$sep_cutoffs <- ests
    return(x)
}