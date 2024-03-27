suppressPackageStartupMessages({
    library(dplyr)
    library(diffcyt)
    library(SummarizedExperiment)
})

set.seed(3004)
data(PBMC_fs, PBMC_panel, PBMC_md)
x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
x <- cluster(x, seed = 3004, verbose = FALSE)

k <- "meta20"
es <- assay(x, "exprs")
kids <- cluster_ids(x, k)

design <- createDesignMatrix(ei(x), cols_design = 2:3)
contrast <- createContrast(c(0, 1, 0, 0, 0))

da <- diffcyt(x,
    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", verbose = FALSE,
    clustering_to_use = k, design = design, contrast = contrast)

ds <- diffcyt(x,
    analysis_type = "DS", method_DS = "diffcyt-DS-limma", verbose = FALSE,
    clustering_to_use = k, design = design, contrast = contrast)

da <- rowData(da$res)
ds <- rowData(ds$res)

test_that("plotDiffHeatmap() - DA", {
    p <- plotDiffHeatmap(x, da, all = TRUE, sort_by = "none", normalize = FALSE)
    expect_is(p, "Heatmap")
    y <- p@matrix
    expect_gte(min(y), 0)
    expect_lte(max(y), 1)
    expect_equivalent(colSums(y), rep(1, ncol(y)))
    expect_equal(dim(y), c(nlevels(kids), nlevels(x$sample_id)))
    expect_identical(rownames(y), levels(kids))
    expect_identical(colnames(y), levels(x$sample_id))
    expect_identical(c(prop.table(table(kids, x$sample_id), 2)), c(y))
})

test_that("plotDiffHeatmap() - DS", {
    p <- plotDiffHeatmap(x, ds, top_n = (n <- 5), 
        sort_by = "padj", lfc = 0, normalize = FALSE)
    expect_is(p, "Heatmap")
    df <- mutate_if(data.frame(ds), is.factor, as.character)
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

test_that("plotDiffHeatmap() - DA; cluster filtering", {
    for (ks in lapply(c(1, 5, 10), sample, x = levels(kids))) {
        y <- filterSCE(x, !cluster_id %in% ks, k = k)
        p <- plotDiffHeatmap(y, da, all = TRUE, sort_by = "none")
        expect_is(p, "Heatmap")
        expect_identical(rownames(p@matrix), setdiff(levels(kids), ks))
    }
})

test_that("plotDiffHeatmap() - DS; cluster filtering", {
    for (ks in lapply(c(1, 3, 5), sample, x = levels(kids))) {
        y <- filterSCE(x, !cluster_id %in% ks, k = k)
        nk <- nlevels(cluster_ids(y, k))
        p <- plotDiffHeatmap(y, ds, all = TRUE, top_n = Inf)
        expect_is(p, "Heatmap")
        hm_ks <- gsub(".*\\(([0-9]+)\\)", "\\1", rownames(p@matrix))
        expect_true(!any(ks %in% hm_ks))
    }
})

test_that("plotDiffHeatmap() - DS; marker filtering", {
    for (ms in lapply(c(1, 3, 5), sample, x = state_markers(x))) {
        n <- sample(nlevels(kids)-5, 1)
        ks <- sample(levels(kids), n)
        y <- filterSCE(x[ms, ], cluster_id %in% ks, k = k)
        p <- plotDiffHeatmap(y, ds, all = TRUE, top_n = Inf)
        expect_is(p, "Heatmap")
        hm_ms <- gsub("(.*)\\(.*", "\\1", rownames(p@matrix))
        expect_true(all(hm_ms %in% ms))
    }
})
