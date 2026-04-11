library(SingleCellExperiment)
data(PBMC_fs, PBMC_panel, PBMC_md)
x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)

test_that("runDR()", {
    expect_error(runDR(x, "x"))
    x <- runDR(x, "PCA", cells=(n <- 20))
    dr <- reducedDim(x, "PCA")
    expect_is(x, "SingleCellExperiment")
    expect_true(reducedDimNames(x) == "PCA")
    expect_identical(sum(!is.na(dr))/ncol(dr), n*nlevels(x$sample_id))
    expect_silent(runDR(x[, !is.na(dr[, 1])], "UMAP", reddim.type="PCA"))
})

test_that("runDR() - use specific subset of features", {
    i <- sample(rownames(x), 10)
    cs <- sample(seq_len(ncol(x)), 100)
    dr1 <- runTsne(assay(x[i, cs], "exprs"))
    dr2 <- reducedDim(runDR(x[, cs], "TSNE", features=i))
    expect_identical(dr1, dr2)
})

test_that("runDR() - use type/state features only", {
    for (c in c("type", "state")) {
        i <- rownames(x)[rowData(x)$marker_class == c]
        set.seed(1); dr1 <- reducedDim(runDR(x, "TSNE", features=i, cells=50))
        set.seed(1); dr2 <- reducedDim(runDR(x, "TSNE", features=c, cells=50))
        # check that same subset of cells has been used
        expect_identical(
            (cs1 <- which(!is.na(dr1[, 1]))), 
            (cs2 <- which(!is.na(dr2[, 1]))))
        # run using 'scrapper' as reference
        cs <- sample(seq_len(ncol(x)), 100)
        expect_equal(
            runTsne(assay(x[i, cs], "exprs")),
            reducedDim(runDR(x[i, cs], "TSNE", features=NULL)))
        dr <- reducedDim(runPca.se(x[, cs], features=i, assay.type="exprs"))
        set.seed(1); expect_equal(reducedDim(runDR(x[i, cs], "PCA", features=NULL)), dr)
        set.seed(1); expect_equal(reducedDim(runDR(x[, cs], "PCA", features=i)), dr)
    }
})
