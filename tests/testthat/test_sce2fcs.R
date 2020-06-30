library(flowCore)
library(SingleCellExperiment)

data(sample_ff, sample_key)
data(PBMC_fs, PBMC_panel, PBMC_md)
ids <- rownames(sample_key)

test_that("sce2fcs()", {
    # split by barcode population
    x <- prepData(sample_ff, by_time = FALSE)
    x <- assignPrelim(x, sample_key, verbose = FALSE)
    expect_error(sce2fcs(x, split_by = "x"))
    expect_error(sce2fcs(x, split_by = "bc_id", assay = "x"))
    expect_is(sce2fcs(x, split_by = NULL), "flowFrame")
    a <- sample(assayNames(x), 1)
    fs <- sce2fcs(x, assay = a, split_by = "bc_id")
    expect_is(fs, "flowSet")
    expect_equivalent(fsApply(fs, nrow), c(table(x$bc_id)))
    m <- grep(id <- sample(ids, 1), fsApply(fs, identifier))
    expect_equivalent(tolerance = 1e-6,
        t(exprs(fs[[m]])), assay(x, a)[, x$bc_id == id])
    # missing populations
    ids_ex <- sample(ids, 5)
    ids_in <- setdiff(ids, ids_ex)
    x <- x[, !x$bc_id %in% ids_ex]
    y <- sce2fcs(x, split_by = "bc_id")
    expect_equivalent(gsub(".*\\.", "", fsApply(y, identifier)), ids_in)
    # split by cluster assignment
    x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
    x <- cluster(x, maxK = 5, verbose = FALSE)
    x$cluster_id <- cluster_ids(x, "meta5")
    y <- sce2fcs(x, split_by = "cluster_id")
    expect_equivalent(fsApply(y, nrow), c(table(x$cluster_id)))
    # with propagation of dimension reductions
    x <- runDR(x, dr = "PCA", ncomponents = 2)
    y <- sce2fcs(x, split_by = "sample_id", keep_dr = TRUE)
    expect_true(length(y) == nlevels(x$sample_id))
    expect_true(ncol(y[[1]]) == nrow(x) + ncol(reducedDim(x)))
})
