suppressPackageStartupMessages({
    library(flowCore)
    library(SingleCellExperiment)
})

data(GvHD)
fs <- GvHD[seq_len(5)]
sce <- prepData(fs[[1]])

test_that(".read_fs() - input arguments", {
    u <- fs
    class(u) <- "x"
    expect_error(.read_fs(u))
    u <- flowSet_to_list(fs)
    u[[1]] <- "x"
    expect_error(.read_fs(u))
    u <- system.file("extdata", package = "flowCore")
    u <- file.path(u, dir(u))[seq_len(3)]
    u[1] <- "x"
    expect_error(.read_fs(u))
})
test_that(".read_fs() - 'flowSet' as input", {
    x <- .read_fs(fs)[[1]]
    expect_identical(x, fs)
})
test_that(".read_fs() - 'flowFrame's as input", {
    ffs <- list(fs[[1]], fs[[2]])
    x <- .read_fs(ffs)[[1]]
    expect_equal(x, flowSet(ffs))
})
test_that(".read_fs() - FCS files as input", {
    fcs <- system.file("extdata", package = "flowCore")
    fcs <- file.path(fcs, dir(fcs))[seq_len(3)]
    x <- .read_fs(fcs)[[1]]
    y <- lapply(fcs, read.FCS,
        transformation = FALSE,
        truncate_max_range = FALSE)
    expect_is(x, "flowSet")
    expect_equal(x, flowSet(y))
})

# ------------------------------------------------------------------------------

test_that(".transform() - input arguments", {
    cfs <- seq_len(nrow(sce))
    names(cfs) <- channels(sce)
    cfs1 <- unname(cfs)
    cfs2 <- cfs
    names(cfs2)[1] <- "x"
    for (cf in list(0, "x", cfs1, cfs2)) 
        expect_error(.transform(sce, cf))
})
test_that(".transform() - single cofactor", {
    cf <- sample(10, 1)
    cfs <- rep(cf, nrow(sce))
    names(cfs) <- channels(sce)
    x <- .transform(sce, cf)
    y <- .transform(sce, cfs)
    int_metadata(x)$cofactor <- NULL
    int_metadata(y)$cofactor <- NULL
    expect_identical(x, y)
})
test_that(".transform() - channel-specific cofactors", {
    cfs <- seq_len(nrow(sce))
    names(cfs) <- channels(sce)
    expect_identical(
        .transform(sce, cfs),
        .transform(sce, sample(cfs)))
})
