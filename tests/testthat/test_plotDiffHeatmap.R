suppressPackageStartupMessages({
    library(purrr)
    library(dplyr)
    library(diffcyt)
    library(data.table)
    library(reshape2)
})

data(PBMC_fs, PBMC_panel, PBMC_md)
x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
x <- cluster(x, verbose = FALSE)
es <- assay(x, "exprs")
kids <- cluster_ids(x, (k <- "meta20"))

design <- createDesignMatrix(PBMC_md, cols_design = 3:4)
contrast <- createContrast(c(0, 1, 0, 0, 0))

da <- diffcyt(x, clustering_to_use = k, design = design, contrast = contrast,
    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", verbose = FALSE)

ds <- diffcyt(x, clustering_to_use = k, design = design, contrast = contrast,
    analysis_type = "DS", method_DS = "diffcyt-DS-limma", verbose = FALSE)

test_that("plotDiffHeatmap() - DA", {
    p <- plotDiffHeatmap(x, da, order = FALSE, all = TRUE, normalize = FALSE)
    expect_is(p, "Heatmap")
    y <- p@matrix
    expect_gte(min(y), 0)
    expect_lte(max(y), 1)
    expect_true(all(colSums(y) == 1))
    expect_equal(dim(y), c(nlevels(kids), nlevels(x$sample_id)))
    expect_identical(rownames(y), levels(kids))
    expect_identical(colnames(y), levels(x$sample_id))
    expect_identical(c(prop.table(table(kids, x$sample_id), 2)), c(y))
})

test_that("plotDiffHeatmap() - DS", {
    p <- plotDiffHeatmap(x, ds, top_n = (n <- 10),
        order = TRUE, row_anno = FALSE, normalize = FALSE)
    expect_is(p, "Heatmap")
    df <- data.frame(rowData(ds$res))
    df <- mutate_if(df, is.factor, as.character)
    df <- df[order(df$p_adj)[seq_len(n)], ]
    y <- p@matrix
    expect_equal(dim(y), c(n, nlevels(x$sample_id)))
    expect_identical(colnames(y), levels(x$sample_id))
    expect_identical(gsub("\\(.*", "", rownames(y)), df$marker_id)
    replicate(10, {
        i <- sample(nrow(df), 1)
        s <- sample(levels(x$sample_id), 1)
        j <- sprintf("%s(%s)", df$marker_id, df$cluster_id)[i]
        cs <- x$sample_id == s & kids == df$cluster_id[i]
        expect_identical(median(es[df$marker_id[i], cs]), p@matrix[j, s])
    })
})

test_that("plotMultiHeatmap() - DA; filtering", {
    for (ks in lapply(c(1, 5, 10), sample, x = levels(kids))) {
        y <- filterSCE(x, !cluster_id %in% ks, k = k)
        p <- plotDiffHeatmap(y, da, all = TRUE, order = FALSE)
        expect_is(p, "Heatmap")
        expect_identical(rownames(p@matrix), setdiff(levels(kids), ks))
    }
})

test_that("plotMultiHeatmap() - DS; filtering", {
    for (ks in lapply(c(1, 3, 5), sample, x = levels(kids))) {
        y <- filterSCE(x, !cluster_id %in% ks, k = k)
        nk <- nlevels(cluster_ids(y, k))
        p <- plotDiffHeatmap(y, ds, top_n = (n <- 2)*nk,
            order = FALSE, row_anno = FALSE)
        expect_is(p, "Heatmap")
        expect_identical(
            rep(setdiff(levels(kids), ks), 2),
            gsub(".*\\(([0-9]+)\\)", "\\1", rownames(p@matrix)))
    }
    for (ms in lapply(c(1, 3, 5), sample, x = state_markers(x))) {
        n <- sample(nlevels(kids)-5, 1)
        ks <- sample(levels(kids), n)
        y <- filterSCE(x[ms, ], cluster_id %in% ks, k = k)
        p <- plotDiffHeatmap(y, ds, all = TRUE, order = FALSE)
        expect_is(p, "Heatmap")
        expect_identical(
            table(rep(ms, n)),
            table(gsub("(.*)\\(.*", "\\1", rownames(p@matrix))))
    }
})
