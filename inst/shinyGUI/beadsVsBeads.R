library(flowCore)
library(CATALYST)
ff <- concatFCS(raw_data)
ms <- gsub("[[:punct:][:alpha:]]", "", colnames(ff))
beadCols <- get_bead_cols(colnames(ff), "dvs")
beadMs <- ms[beadCols]
n_beads <- length(beadCols)
key <- data.frame(matrix(c(0, 0, rep(1, n_beads)), ncol=2+n_beads,
                         dimnames=list("is_bead", c(191, 193, beadMs))), check.names=FALSE)
beadInds <- get_bead_inds(ff, key)
es <- exprs(ff)[, beadCols]
data <- es[beadInds, ]
mhlDists <- sqrt(stats::mahalanobis(
    x=es, center=colMeans(data), cov=stats::cov(data)))

mu <- colMeans(exprs(ff)[beadInds, beadCols])
sigma <- cov(exprs(ff)[beadInds, beadCols])
x <- exprs(ff)[, beadCols]
# (Y(I,:)-mu)*inv(SIGMA)*(Y(I,:)-mu)'
mhlDists <- sqrt(sapply(seq_len(nrow(x)), function(i) (x[i, ] - mu) %*% solve(sigma) %*% (x[i, ] - mu)))
tmp <- t(t(x) - mu) %*% solve(sigma) %*% t(t(t(x) - mu))
tmp <- t(t(x) - mu) %*% solve(sigma) %*% t(t(x) - mu)
tmp <- sapply(seq_len(nrow(tmp)), function(i) tmp[i, upper.tri(tmp)[i, ]])

dists <- tmp[upper.tri(tmp)]

get_axes <- function(df, cf) {
    min <- max(apply(df, 2, function(i) -ceiling(abs(min(i))*2)/2))
    max <- max(apply(df, 2, function(i) ceiling(max(i)*2)/2))
    
    if (min != 0) {
        tcks <- c(-10^(ceiling(log10(abs(sinh(min)*cf))):0), 
            0, 10^(0:ceiling(log10(sinh(max)*cf))))
    } else {
        tcks <- c(0,  10^(0:ceiling(log10(sinh(max)*cf))))
    }
    labs <- parse(text=gsub("[[:digit:]]*e", " 10^",
        format(tcks, scientific=TRUE)))
    labs[tcks == 0] <- ""
    tcks <- asinh(tcks/cf)
    return(list(tcks, labs))
}

beadsVsBeads <- function(beads, mhlDists) {
    
    # subsample for visualization 
    N <- nrow(beads)
    if (N > 1e5) {
        inds <- sample(N, 1e5)
        beads <- beads[inds, ]
        mhlDists <- mhlDists[inds]
    }
    
    chs <- colnames(beads)
    tbeads <- asinh(beads/5)

    n <- ncol(beads)
    if (n %% 2 != 0) {
        combis <- matrix(c(1:n, n-1), nrow=n/2)
    } else {
        combis <- matrix(1:n, ncol=n/2)
    }

    nPlots <- ncol(combis)
    nEvents <- nrow(data)
    maxDist <- ceiling(max(mhlDists))
    if (maxDist < 30) {
        stp <- 1 
    } else if (maxDist < 100) {
        stp <- 5
    } else {
        stp <- 10
    }
    
    # get axis limits and labels
    temp <- get_axes(tbeads, 5)
    tcks <- temp[[1]]
    labs <- temp[[2]]

    first <- TRUE
    ps <- vector("list", nPlots+1)
    for (i in seq_len(nPlots)) {
        df <- data.frame(
            x=tbeads[, combis[1,i]],
            y=tbeads[, combis[2,i]],
            z=mhlDists)
        ps[[i]] <- ggplot(df) + 
            geom_point(size=.5, aes_string(x="x", y="y", color="z")) +
            geom_vline(xintercept=0, size=.5) +
            geom_hline(yintercept=0, size=.5) +
            labs(x=chs[combis[1,i]], y=chs[combis[2,i]]) +
            coord_cartesian(xlim=tcks, ylim=tcks, expand=.1) +
            scale_x_continuous(breaks=tcks, labels=labs) +
            scale_y_continuous(breaks=tcks, labels=labs) +
            scale_color_gradientn(
                colours=c(RColorBrewer::brewer.pal(11, "Spectral")),
                limits=c(0, maxDist), breaks=seq(0, maxDist, stp),
                name="Distance from identified beads") +
            guides(colour=FALSE) + theme_classic() + theme(
                aspect.ratio=1,
                plot.margin=unit(c(0, 0, .25, .25), "cm"),
                panel.grid.minor=element_blank(),
                panel.grid.major=element_blank(),
                axis.ticks=element_line(size=.5),
                axis.title=element_text(size=14, face="bold"),
                axis.text.x=element_text(size=12, color="black", angle=30, hjust=1, vjust=1),
                axis.text.y=element_text(size=12, color="black", angle=30))
        if (first) {
            ps[[i]] <- ps[[i]] + 
                guides(colour=guide_colourbar(
                    title.position="top", title.hjust=.5)) + 
                theme(
                    legend.direction="horizontal",
                    legend.title=element_text(face="bold", size=12),
                    legend.text=element_text(size=10),
                    legend.key=element_rect(colour="white"),
                    legend.key.height=unit(.5, "cm"),
                    legend.key.width=unit(.55/nPlots, "npc"))
            lgd <- get_legend(ps[[i]])
            ps[[i]] <- ps[[i]] + guides(colour=FALSE)
            first <- FALSE
        }
    }
    ps[[i+1]] <- lgd
    grid.arrange(grobs=ps, 
        layout_matrix=rbind(rep(nPlots+1, nPlots), seq_len(nPlots)), 
        heights=c(2, 8),
        widths=rep(8, nPlots))
}
pdf("/Users/HLC/Desktop/beadsVsBeads2.pdf", 15, 6)
beadsVsBeads(es, mhlDists)
dev.off()
# mhlDists[mhlDists>40] <- 30
