library(scater)
library(SingleCellExperiment)

data(PBMC_fs, PBMC_panel, PBMC_md)
x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)

test_that("runDR()", {
    expect_error(runDR(x, "x"))
    x <- runDR(x, "MDS", cells = (n <- 20))
    dr <- reducedDim(x, "MDS")
    expect_is(x, "SingleCellExperiment")
    expect_true(reducedDimNames(x) == "MDS")
    expect_identical(sum(!is.na(dr))/2, n*nlevels(x$sample_id))
    expect_silent(runDR(x[, !is.na(dr[, 1])], "TSNE", dimred = "MDS"))
})

test_that("runDR() - use type/state features only", {
    for (c in c("type", "state")) {
        i <- rownames(x)[rowData(x)$marker_class == c]
        set.seed(1); dr1 <- reducedDim(runDR(x, "MDS", features = i, cells = 50))
        set.seed(1); dr2 <- reducedDim(runDR(x, "MDS", features = c, cells = 50))
        # check that same subset of cells has been used
        expect_identical(
            (cs1 <- which(!is.na(dr1[, 1]))), 
            (cs2 <- which(!is.na(dr2[, 1]))))
        # run using 'scater' as reference
        cs <- sample(seq_len(ncol(x)), 100)
        expect_equal(
            reducedDim(runDR(x[i, cs], "MDS", features = NULL)),
            reducedDim(runMDS(x[, cs], subset_row = i, exprs_values = "exprs")))
        dr <- reducedDim(runMDS(x[, cs], subset_row = i, exprs_values = "exprs"))
        set.seed(1); expect_equal(reducedDim(runDR(x[i, cs], "MDS", features = NULL)), dr)
        set.seed(1); expect_equal(reducedDim(runDR(x[, cs], "MDS", features = i)), dr)
    }
})

test_that("runDR() - use specific subset of features", {
    i <- sample(rownames(x), 10)
    cs <- sample(seq_len(ncol(x)), 100)
    dr1 <- reducedDim(runDR(x[, cs], "MDS", features = i))
    dr2 <- reducedDim(runMDS(x[, cs], subset_row = i, exprs_values = "exprs"))
    expect_identical(dr1, dr2)
})
