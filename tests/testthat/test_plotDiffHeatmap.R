context("differential")

library(purrr)
library(dplyr)
library(diffcyt)
library(data.table)
library(reshape2)
library(SingleCellExperiment)

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
    foo <- capture.output({
        pdf(NULL)
        p <- plotDiffHeatmap(x, da, top_n = (n <- 10), normalize = FALSE)
        dev.off()
    })
    expect_is(p, "HeatmapList")
    y <- p@ht_list[[2]]@matrix
    o <- order(rowData(da$res)$p_adj)[seq_len(n)]
    expect_gte(min(y), 0); expect_lte(max(y), 1)
    expect_equal(dim(y), c(n, nlevels(x$sample_id)))
    expect_identical(rownames(y), levels(kids)[o])
    expect_identical(colnames(y), levels(x$sample_id))
    expect_identical(c(prop.table(table(kids, x$sample_id), 2)[o, ]), c(y))
})

test_that("plotDiffHeatmap() - DS results", {
    foo <- capture.output({
        pdf(NULL)
        p <- plotDiffHeatmap(x, ds, top_n = (n <- 10), normalize = FALSE)
        dev.off()
    })
    expect_is(p, "HeatmapList")
    df <- data.frame(rowData(ds$res))
    df <- mutate_if(df, is.factor, as.character)
    y <- p@ht_list[[2]]@matrix
    o <- order(df$p_adj)[seq_len(n)]
    expect_equal(dim(y), c(n, nlevels(x$sample_id)))
    expect_identical(rownames(y), as.character(rowData(ds$res)$marker_id[o]))
    expect_identical(colnames(y), levels(x$sample_id))
    cs <- data.table(i = seq_len(ncol(x)), 
        data.frame(colData(x))) %>% 
        do({.$cluster_id = kids; .}) %>% 
        split(by = c("cluster_id", "sample_id"), 
            sorted = TRUE, flatten = FALSE) %>% 
        map_depth(2, "i")
    ms <- apply(df[o, ], 1, function(u) {
        k <- as.character(u[["cluster_id"]])
        m <- as.character(u[["marker_id"]])
        sapply(cs[[k]], function(i) median(assay(x)[m, i]))
    })
    expect_identical(c(t(ms)), c(y))
})
