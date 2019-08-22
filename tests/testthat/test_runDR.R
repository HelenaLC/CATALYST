context("dimension reduction")

suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

data(PBMC_fs, PBMC_panel, PBMC_md)
x <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)

test_that("row subsetting", {
    i <- sample(nrow(x), 100)
    expect_identical(
        as.data.frame.matrix(rowData(x[i, ])),
        as.data.frame.matrix(rowData(x)[i, ]))
})

test_that(".rotate_daf", {
    y <- .rotate_daf(x)
    expect_is(y, "SingleCellExperiment")
    expect_identical(colData(y), rowData(x))
    expect_identical(rowData(y), colData(x))
    expect_identical(metadata(y), metadata(x))
    expect_identical(reducedDims(y), reducedDims(x))
})

test_that("throw error if DR already exists", {
    i <- sample(nrow(x), 10)
    y <- runDR(x, "PCA", i)
    expect_error(runDR(y, "PCA", i, overwrite = FALSE))
    expect_silent(runDR(y, "PCA", i, overwrite = TRUE))
})

test_that("throw warning when 'use_dimred' doesn't exist", {
    i <- sample(nrow(x), 10)
    y <- runDR(x, "PCA", i)
    expect_error(runDR(y, "TSNE", use_dimred = "x"))
    expect_silent(runDR(y, "TSNE", i, use_dimred = "PCA"))
})

test_that("throw error when 'use_dimred' is supplied to PCA", {
    i <- sample(nrow(x), 10)
    y <- runDR(x, "MDS", i)
    expect_error(runDR(y, "PCA", i, use_dimred = "MDS"))
})

test_that("multiple DRs w/ different #rows work", {
    n1 <- 20; n2 <- 50
    i1 <- sample(nrow(x), n1)
    i2 <- sample(nrow(x), n2)
    y <- runDR(x, "PCA", i1)
    y <- runDR(y, "MDS", i2)
    expect_equal(nrow(reducedDim(y, "PCA")), n1)
    expect_equal(nrow(reducedDim(y, "MDS")), n2)
})

test_that("TSNE using PCA from different cells fails", {
    i1 <- seq_len(10)[-1]
    i2 <- seq_len(10)[-2]
    y <- runDR(x, "PCA", i1)
    expect_error(runDR(y, "TSNE", i2, use_dimred = "PCA"))
})

test_that("TSNE using PCA from identical cells succeeds", {
    i <- sample(nrow(x), 10)
    y <- runDR(x, "PCA", i)
    expect_silent(runDR(y, "TSNE", i, use_dimred = "PCA"))
})

test_that("DR on subset i & rows_to_use = i is identical", {
    i <- sample(nrow(x), 10)
    expect_identical(
        reducedDim(runDR(x, "PCA", i)),
        reducedDim(runDR(x[sort(i), ], "PCA", NULL)))
})
