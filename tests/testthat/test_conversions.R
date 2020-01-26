library(flowCore)
library(SingleCellExperiment)

data(sample_ff, sample_key)
data(PBMC_fs, PBMC_panel, PBMC_md)
ids <- rownames(sample_key)

test_that("sce2fcs()", {
    # split by barcode population
    x <- fcs2sce(sample_ff, by_time = FALSE)
    x <- assignPrelim(x, sample_key, verbose = FALSE)
    expect_error(sce2fcs(x, split_by = "x"))
    expect_error(sce2fcs(x, split_by = "bc_id", assay = "x"))
    expect_is(sce2fcs(x, split_by = NULL), "flowFrame")
    expect_is(fs <- sce2fcs(x, split_by = "bc_id"), "flowSet")
    expect_equivalent(fsApply(fs, nrow), c(table(x$bc_id)))
    m <- grep(id <- sample(ids, 1), fsApply(fs, identifier))
    expect_identical(t(exprs(fs[[m]])), assay(x, "exprs")[, x$bc_id == id])
    # missing populations
    ids_ex <- sample(ids, 5)
    ids_in <- setdiff(ids, ids_ex)
    x <- x[, !x$bc_id %in% ids_ex]
    y <- sce2fcs(x, split_by = "bc_id")
    expect_equivalent(gsub(".*\\.", "", fsApply(y, identifier)), ids_in)
    # split by cluster assigment
    x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
    x <- cluster(x, maxK = 5, verbose = FALSE)
    y <- sce2fcs(x, split_by = "cluster_id")
    expect_equivalent(fsApply(y, nrow), c(table(x$cluster_id)))
    # with propagation of dimension reductions
    x <- runDR(x, dr = "MDS", cells = 100)
    y <- sce2fcs(x, split_by = "sample_id", keep_dr = TRUE)
    expect_true(length(y) == nlevels(x$sample_id))
    expect_true(ncol(y[[1]]) == nrow(x) + ncol(reducedDim(x)))
})

test_that("fcs2sce()", {
    # construct artificial flowFrame
    n <- 100; m <- 50
    x <- matrix(runif(n * m), n, m)
    colnames(x) <- seq_len(m)
    ff <- flowFrame(x)
    # dimensions should be reversed
    y <- fcs2sce(ff, by_time = FALSE)
    expect_is(y, "SingleCellExperiment")
    expect_identical(rev(dim(y)), dim(x))
    # data should be unchanged
    expect_identical(t(assay(y)), x)
    # with transformation
    expect_error(fcs2sce(ff, transform = "x"))
    y <- fcs2sce(ff, by_time = FALSE,
        transform = TRUE, cofactor = (cf <- 20))
    expect_true("exprs" %in% assayNames(y))
    expect_identical(assay(y, "exprs"), asinh(t(exprs(ff))/cf))
    # construct artifical flowSet
    i <- sample(seq_len(n), 10)
    j <- seq_len(n)[-i]
    fs <- flowSet(ff[i, ], ff[j, ])
    # should render message when 
    # acquisition times aren't available
    expect_message(fcs2sce(fs))
    # assure flowFrames are ordered correctly
    description(fs[[1]])$`$BTIM` <- 2
    description(fs[[2]])$`$BTIM` <- 1
    expect_silent(y <- fcs2sce(fs))
    expect_equal(assay(y)[, seq_along(j)], t(exprs(ff[j, ])))
    expect_equal(assay(y)[, seq_along(i)+length(j)], t(exprs(ff[i, ])))
})
