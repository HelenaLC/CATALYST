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
kids <- cluster_ids(x, (k <- "meta20"))

design <- createDesignMatrix(PBMC_md, cols_design = 3:4)
contrast <- createContrast(c(0, 1, 0, 0, 0))

da <- diffcyt(x, clustering_to_use = k, design = design, contrast = contrast, 
    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", verbose = FALSE)

ds <- diffcyt(x, clustering_to_use = k, design = design, contrast = contrast, 
    analysis_type = "DS", method_DS = "diffcyt-DS-limma", verbose = FALSE)

test_that("plotDiffHeatmap() - DA results", {
    p <- plotDiffHeatmap(x, da, order = FALSE, all = TRUE, normalize = FALSE)
    expect_is(p, "HeatmapList")
    y <- p@ht_list[[2]]@matrix
    expect_gte(min(y), 0)
    expect_lte(max(y), 1)
    expect_equivalent(rowSums(y), rep(1, nrow(y)))
    expect_equal(dim(y), c(nlevels(kids), nlevels(x$sample_id)))
    expect_identical(rownames(y), levels(kids))
    expect_identical(colnames(y), levels(x$sample_id))
    expect_identical(c(prop.table(table(kids, x$sample_id), 1)), c(y))
})

test_that("plotDiffHeatmap() - DS results", {
    p <- plotDiffHeatmap(x, ds, top_n = (n <- 10), 
        order = TRUE, row_anno = FALSE)
    expect_is(p, "HeatmapList")
    df <- data.frame(rowData(ds$res))
    df <- mutate_if(df, is.factor, as.character)
    df <- df[order(df$p_adj)[seq_len(n)], ]
    y <- p@ht_list[[2]]@matrix
    expect_equal(dim(y), c(n, nlevels(x$sample_id)))
    expect_identical(rownames(y), df$marker_id)
    expect_identical(colnames(y), levels(x$sample_id))
})
