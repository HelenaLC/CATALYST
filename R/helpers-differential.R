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
