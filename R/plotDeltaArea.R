# ==============================================================================
# Relative change in area under CDF curve
# ------------------------------------------------------------------------------

#' @rdname plotDeltaArea
#' @title Relative change in are under CDF curve
#'
#' @param x a list as returned by \code{\link{ConsensusClusterPlus}}.
#' 
#' @return a \code{\link{ggplot}} object.
#' 
#' @details 
#' The delta area represents the amount of extra cluster stability gained when 
#' clustering into k groups as compared to k-1 groups. It can be expected that 
#' high stability of clusters can be reached when clustering into the number of 
#' groups that best fits the data. The "natural" number of clusters present in 
#' the data should thus corresponds to the value of k where there is no longer 
#' a considerable increase in stability (pleateau onset).
#' 
#' @examples
#' 
#' @author
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' @references 
#' Nowicka M, Krieg C, Weber LM et al. 
#' CyTOF workflow: Differential discovery in 
#' high-throughput high-dimensional cytometry datasets.
#' \emph{F1000Research} 2017, 6:748 (doi: 10.12688/f1000research.11622.1)
#' 
#' @import FlowSOM ggplot2
#' @importFrom flowCore flowSet fsApply
#' @importFrom graphics hist
#' @importFrom matrixStats colQuantiles
#' @importFrom reshape2 melt
#' @export
# ==============================================================================

setMethod(f="plotDeltaArea", 
    signature=signature(x="list"), 
    definition=function(x) {
        # check that input corresponds to list
        # as returned by ConsensusClusterPlus 
        is_ConsensusClusterPlus_list(x)
        
        k <- length(x)
        areaK <- numeric(k-1)
        v <- lapply(mc[seq_len(k)[-1]], function(x) triangle(x$ml))
        # empirical CDF distribution
        h <- lapply(v, function(x) {
            h <- graphics::hist(x, breaks=seq(0, 1, .01), plot=FALSE)
            h$counts <- cumsum(h$counts) / sum(h$counts)
            return(h)
        })
        # calculate area under CDF curve, by histogram method
        areaK <- sapply(h, function(x) cumsum(x$counts * .01)[100])
        
        # calculate proportional increase relative to prior k
        # w/ initial AUC at k = 2
        deltaK <- c(areaK[1], diff(areaK) / areaK[seq_len(k-2)])
        
        # plot relative change in area under CDF curve vs. k
        df <- data.frame(k=1+seq_len(k-1), y=deltaK)
        max <- ceiling(max(df$y)*2)/2
        ggplot(df, aes(x=k, y=y)) + geom_line(color="steelblue") + 
            geom_point(size=2.5, color="navy") + coord_fixed(4) +
            scale_x_continuous(breaks=seq(2, 20, 2), expand=c(0,.5)) +
            scale_y_continuous(limits=c(0, max), expand=c(0,.125), 
                breaks=function(x) seq(x[1]+.125, x[2], .5)) +
            ylab("Relative change in area under CDF curve") +
            ggtitle("Delta area") + theme_classic() + 
            theme(plot.title=element_text(face="bold"),
                axis.text=element_text(color="black"),
                panel.grid.major=element_line(color="lightgrey", size=.25))
    }
)