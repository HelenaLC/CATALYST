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
#' @param verbose 
#' logical. Should extra information on progress be reported? Defaults to TRUE.
#'
#' @return
#' Will update the \code{sep_cutoffs}, \code{mhl_cutoff}, \code{counts} and 
#' \code{yields} slots of the input \code{\link{dbFrame}} and return the latter.
#' 
#' @examples
#' data(sample_ff, sample_key)
#' re <- assignPrelim(x = sample_ff, y = sample_key)
#' estCutoffs(x = re)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @importFrom stats loess
#' @importFrom drc drm LL.4

# ------------------------------------------------------------------------------

setMethod(f="estCutoffs", 
    signature=signature(x="dbFrame"), 
    definition=function(x, verbose=TRUE) {
        
        cutoffs <- seq(0, 1, .01)
        n_bcs <- nrow(bc_key(x))
        ests <- numeric(n_bcs)
        
        # f0 <- function(x, b, c, d, e) c+(d-c)/(1+exp(b*(log(x)-log(e))))
        # f1 <- function(x, b, c, d, e) b*e^b*x^(b-1)*(c-d)/(x^b+e^b)^2
        # f2 <- function(x, b, c, d, e) b*e^b*x^(b-2)*((b-1)*e^b-(b+1)*x^b)*(c-d)/(x^b+e^b)^3
        # f3 <- function(x, b, c, d, e) (b*e^b*x^(b-3)*(b^2+3*b+2)*x^(2*b)-4*(b^2-1)*e^b*x^b+(b^2-3*b+2)*e^(2*b)*(c-d))/(x^b+e^b)^4
        f4 <- function(x, b, c, d, e) (b*(c-d)*e^b*x^(-4+b)*((-6+11*b-6*b^2+b^3)*e^(3*b)+(-18+11*b+18*b^2-11*b^3)*e^(2*b)*x^b+(-18-11*b+18*b^2+11*b^3)*e^b*x^(2*b)-(6+11*b+6*b^2+b^3)*x^(3*b)))/(e^b+x^b)^5
        # root <- function(b, c, d, e) 
        #     ((2*e^b*b^2-sqrt(3)*sqrt(b^4*e^(2*b)-b^2*e^(2*b))
        #         -2*e^b)/(b^2+3*b+2))^(1/b)
        
        for (i in seq_len(n_bcs)) {
            df <- data.frame(cutoff=cutoffs, yield=as.vector(yields(x)[i, ]))
            fit <- drc::drm(yield~cutoff, data=df, fct=drc::LL.4(fixed=c(NA, 0, NA, NA)))
            #ests[i] <- round(root(fit$coefficients[1], 0, 
            #    fit$coefficients[2], fit$coefficients[3])*100)/100
            b <- fit$coefficients[1]
            c <- 0#fit$coefficients[2]
            d <- fit$coefficients[2]
            e <- fit$coefficients[3]
            root <- which(diff(sign(diff(f4(cutoffs, b, c, d, e)))) == 2) + 1
            ests[i] <- cutoffs[root[1]]
            #ests[i] <- round(root(fit$coefficients[1], fit$coefficients[3])*100)/100
            # ests[i] <- round(root(
            #     fit$coefficients[1], fit$coefficients[2], 
            #     fit$coefficients[3], fit$coefficients[4])*100)/100
        }
        
        names(ests) <- rownames(bc_key(x))
        sep_cutoffs(x) <- ests
        x
    })