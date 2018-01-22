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
    qs <- matrixStats::colQuantiles(x, probs=c(.01, .99))
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
