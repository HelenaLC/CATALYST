# ==============================================================================
# Estimate separation cutoffs
# ------------------------------------------------------------------------------

#' @rdname estTrim
#' @title Estimation of optimal trim value
#' 
#' @description
#' Estimates a trim value that will minimize the sum over squared
#' popultion- and channel-wise squared medians upon compensation.
#'
#' @param x        a \code{\link{dbFrame}}.
#' @param min,max,step 
#'                 specifies sequence of trim values for which 
#'                 compensation should be evaluated.
#' @param out_path specifies in which location output plot is to be generated. 
#'                 Defaults to NULL.
#' @param name_ext a character string. If specified, will be appended 
#'                 to the output plot's name. Defaults to NULL.
#' 
#' @return
#' For each value along \code{seq(min, max, step)}, \code{estTrim} will call
#' \code{\link{computeSpillmat}} with \code{method = "mean"} and the respective 
#' trim parameter. Returned will be an interactive plot displaying channel-wise 
#' medians upon compensation, and the mean squared deviation from 0. Each point 
#' is labeled with the respective interacting channels.
#' 
#' @examples
#' # get single-stained control samples
#' data(ss_exp)
#' 
#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#' 
#' re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' estTrim(x = re, min = 0.02, max = 0.14, step = 0.02)
#'
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @import ggplot2
#' @importFrom graphics plot
#' @importFrom plotly ggplotly
#' @export

# ------------------------------------------------------------------------------

setMethod(f="estTrim", 
    signature=signature(x="dbFrame"), 
    definition=function(x, min = 0.05, max = 0.20, step = 0.01, 
        out_path = NULL, name_ext = NULL) {
        
        trms <- seq(min, max, step)
        ids <- as.numeric(rownames(bc_key(x)))
        nms <- paste(colnames(exprs(x)))
        ms <- gsub("[[:alpha:][:punct:]]", "", nms)
        bc_cols <- which(ms %in% ids)
        ms <- ms[bc_cols]
        mets <- gsub("[[:digit:]]+Di", "", nms[bc_cols])
        spill_cols <- get_spill_cols(as.numeric(ms), mets)
        x@exprs <- exprs(x)[, bc_cols]

        # compute spillover and compensation matrix,
        # and compensate data for each trim value
        sm <- lapply(trms, function(val) computeSpillmat(x, trim=val))
        sm <- lapply(sm, make_symetric)
        cm <- lapply(sm, solve)
        comped <- lapply(cm, function(mat) exprs(x) %*% mat)

        # compute channel-wise medians and mean squared
        # devation from zero for each compensated data
        medians <- lapply(comped, function(data) {
            unlist(sapply(ids, function(id) {
                cols <- spill_cols[[which(ms == id)]]
                if (length(cols) == 1) {
                    median(data[bc_ids(x) == id, cols])
                } else {
                    matrixStats::colMedians(data[
                        bc_ids(x) == id, cols])
                }
            }))
        })
        mse <- vapply(medians, function(m) mean(m^2), numeric(1))

        nTrms <- length(trms)
        ns <- lapply(spill_cols, length)
        n <- sum(unlist(ns))
        spillers <- sapply(seq_along(ns), 
            function(i) rep(nms[bc_cols][i], ns[i]))
        receivers <- lapply(spill_cols, function(i) nms[bc_cols][i])
        
        df <- data.frame(
            Spiller=rep(unlist(spillers), nTrms),
            Receiver=rep(unlist(receivers), nTrms), 
            m=unlist(medians), 
            t=rep(trms, each=n), 
            e=rep(mse, each=n))
        
        xMin <- trms[1]-step
        xMax <- trms[nTrms]+step
        yMin <- floor(  min(df$m)/.5)*.5
        yMax <- ceiling(max(df$m)/.5)*.5
        rect <- data.frame(x1=xMin, x2=xMax, y1=yMax+.4, y2=yMax+.6)
        text <- data.frame(x=trms, y=yMax+.5, e=round(mse, 4))
            
        p <- plot_estTrim(df, trms, xMin, xMax, yMin, yMax, rect, text)
        if (!is.null(out_path)) {
            ggsave(file.path(out_path, paste0("estTrim", name_ext, ".pdf")), 
                plot=p, width=nTrms, height=8)
        } else {
            ggplotly(p, tooltip=c("group", "fill"))
        }
    })
