context("daFrame-subsetting")

library(SingleCellExperiment)

data(PBMC_fs, PBMC_panel, PBMC_md)
x <- daFrame(PBMC_fs, PBMC_panel, PBMC_md)

test_that("using nothing", {
    expect_equal(x[ ], x)
    expect_equal(x[, ], x)
    expect_error(x[, , ])
})

test_that("using logicals", {
    expect_equal(dim(x[TRUE, ]), dim(x))
    expect_equal(dim(x[, TRUE]), dim(x))
    expect_equal(dim(x[TRUE, TRUE]), dim(x))
    
    expect_equal(dim(x[FALSE, ]), c(0, ncol(x)))
    expect_equal(dim(x[, FALSE]), c(nrow(x), 0))
    expect_equal(dim(x[FALSE, FALSE]), c(0, 0))
    
    expect_equal(dim(x[TRUE, FALSE]), c(nrow(x), 0))
    expect_equal(dim(x[FALSE, TRUE]), c(0, ncol(x)))

    l <- as.logical(c(0, 1))
    ni <- sample(seq(2, nrow(x) - 1), 1)
    nj <- sample(seq(2, ncol(x) - 1), 1)
    i <- sample(l, ni, TRUE)
    j <- sample(l, nj, TRUE)
    expect_error(x[i, ])
    expect_error(x[, j])
    expect_error(x[i, j])
})

test_that("using numerics", {
    ni <- 100
    nj <- 10
    i <- sample(nrow(x), ni)
    j <- sample(ncol(x), nj)
    
    expect_equal(nrow(x[i, ]), ni)
    expect_equal(ncol(x[, j]), nj)
    expect_equal(dim(x[i, j]), c(ni, nj))
    
    expect_identical(rowData(x[, j]), rowData(x))
    expect_identical(colData(x[i, ]), colData(x))
    
    expect_equivalent(rowData(x[i, ]), rowData(x)[i, ])
    expect_identical(colData(x[, j]), colData(x)[j, ])
})

test_that("col-subsetting doesn't affect DR", {
    i <- sample(nrow(x), 100)
    j <- sample(ncol(x), 10)
    y <- runDR(x, "PCA", i)
    expect_identical(
        reducedDim(y, "PCA"),
        reducedDim(y[, j], "PCA"))
})

test_that("row-subsetting doesn't affects DR within bonds", {
    i <- seq_len(10)
    y <- runDR(x, "PCA", i)
    expect_identical(
        reducedDim(y, "PCA"),
        reducedDim(y[i, ], "PCA"))
})

test_that("row-subsetting affects DR out of bonds", {
    n <- 10
    i1 <- sample(nrow(x), n)
    i2 <- sample(i1, n-5)
    y <- runDR(x, "PCA", i1)
    m <- match(sort(i1), i2, nomatch = 0)
    expect_equivalent(
        reducedDim(y, "PCA")[m, ],
        reducedDim(y[i2, ], "PCA"))
})

test_that("drop DR out of bonds", {
    n <- 10
    i1 <- sample(nrow(x), n)
    i2 <- seq_len(nrow(x))
    i2 <- setdiff(i2, i1)
    y <- runDR(x, "PCA", i1)
    expect_message(y[i2, ])
})
