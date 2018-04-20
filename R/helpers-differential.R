# ==============================================================================
# colors for 20 ConensusClusterPlus metaclusters
# ------------------------------------------------------------------------------
cluster_cols <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

# ==============================================================================
# scale expression to values b/w 0 and 1 using 
# low (1%) and high (99%) quantiles as boundaries
# ------------------------------------------------------------------------------
scale_exprs <- function(x) {
    qs <- matrixStats::colQuantiles(as.matrix(x), probs=c(.01, .99))
    x_scaled <- t((t(x) - qs[, 1]) / (qs[, 2] - qs[, 1]))
    x_scaled[x_scaled < 0] <- 0
    x_scaled[x_scaled > 1] <- 1
    return(x_scaled)
}

# ------------------------------------------------------------------------------
# validity check for plotNRS: 
# check that input corresponds to list
# as returned by ConsensusClusterPlus 
# ------------------------------------------------------------------------------
is_ConsensusClusterPlus_list <- function(x) {
    slots <- vapply(x[-1], names, character(5))
    valid <- apply(slots, 2, all.equal, 
        c(paste0("consensus", c("Matrix", "Tree", "Class")), "ml", "clrs"))
    valid <- sum(valid) == length(x) - 1
    if (!valid)
        stop("Invalid input: x should be a list 
            as returned by 'ConsensusClusterPlus'.", call.=FALSE)
}

# ==============================================================================
# calculate non-redundancy score (NRS) for ea. sample
# ------------------------------------------------------------------------------
nrs <- function(x, n=3) {
    pc <- prcomp(x, center=TRUE, scale.=FALSE)
    scores <- rowSums(outer(rep(1, ncol(x)),
        pc$sdev[seq_len(n)]^2) * abs(pc$rotation[, seq_len(n)]))
    return(scores)
}

# ==============================================================================
# wrapper for ComplexHeatmap row annotations
# (called by plotClusterHeatmap)
# ------------------------------------------------------------------------------
row_anno <- function(anno, cols, name, clustering, dend) {
    Heatmap(matrix=anno, col=cols, name=name, 
        rect_gp=gpar(col="white"), width=unit(.5, "cm"),
        cluster_rows=clustering, cluster_columns=FALSE,
        show_row_dend=dend, row_dend_reorder=FALSE)
}

# ==============================================================================
# change in area under CDF curve
# ------------------------------------------------------------------------------
triangle <- function(m) {
    n <- ncol(m)
    nm <- matrix(0, ncol=n, nrow=n)
    fm <- m
    nm[upper.tri(nm)] <- m[upper.tri(m)]
    fm <- t(nm) + nm
    diag(fm) <-  diag(m)
    nm <- fm
    nm[upper.tri(nm)] <- NA
    diag(nm) <- NA
    m[lower.tri(nm)]
}

plot_delta_area <- function(mc) {
    # empirical CDF distribution
    maxK <- length(mc)
    v <- lapply(mc[seq_len(maxK)[-1]], function(x) triangle(x$ml))
    h <- lapply(v, function(x) {
        h <- graphics::hist(x, breaks=seq(0, 1, .01), plot=FALSE)
        h$counts <- cumsum(h$counts) / sum(h$counts)
        return(h)
    })
    # calculate area under CDF curve, by histogram method &
    # calculate proportional increase relative to prior k
    areaK <- sapply(h, function(x) cumsum(x$counts * .01)[100])
    deltaK <- c(areaK[1], diff(areaK) / areaK[seq_len(maxK-2)])
    
    df <- data.frame(k=seq_len(maxK)[-1], y=deltaK)
    y_max <- ceiling(max(df$y)*2)/2
    ggplot(df, aes_string(x="k", y="y")) + 
        theme_classic() + geom_line(color="steelblue", lty=2) + 
        geom_point(size=2.5, color="navy") + coord_fixed(4) +
        scale_x_continuous(breaks=seq(2, 20, 2), expand=c(0,.5)) +
        scale_y_continuous(limits=c(0, y_max), expand=c(0,.125), 
            breaks=function(x) seq(x[1]+.125, x[2], .5)) +
        ylab("Relative change in area under CDF curve") +
        theme(plot.title=element_text(face="bold"),
            axis.text=element_text(color="black"),
            panel.grid.major=element_line(color="grey", size=.2))
}

# ==============================================================================
# get differential analysis type from differential test result 
# as returned by 'diffcyt::testDA_*()' & 'diffcyt::testDS_*()'
# ------------------------------------------------------------------------------
get_dt_type <- function(x) {
    k <- length(unique(x$cluster_id))
    da <- c("cluster_id", "logFC", "logCPM", "LR", "p_val", "p_adj")
    ds <- c("cluster_id", "marker", "ID", "logFC", "AveExpr", "t", "p_val", "p_adj", "B")
    if (nrow(x) == k && all.equal(colnames(x), da)) {
        return("DA")
    } else if (nrow(x) == (k * nlevels(x$marker))
        && all.equal(colnames(x), ds)) {
        return("DS")
    } else {
        stop(deparse(substitute(x)), " does not seem to be ", 
            "a valid differential test result.\n",
            "Should be a 'SummarizedExperiment' as returned by ", 
            "'diffcyt::testDA_*()' or 'diffcyt::testDS_*()'.")
    }
}

# ==============================================================================
# wrapper for heatmaps by plotDiffHeatmap()
# ------------------------------------------------------------------------------
diff_hm <- function(matrix, col, name, xlab, ...) {
    Heatmap(matrix, col, name, ..., cluster_columns=FALSE,
        column_title=xlab, column_title_side="bottom",
        clustering_distance_rows="euclidean",
        clustering_method_rows="median",
        column_names_gp=gpar(fontsize=8),
        rect_gp=gpar(col='white'))
}
