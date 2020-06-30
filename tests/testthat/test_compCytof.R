data(ss_exp)
library(flowCore)
library(SingleCellExperiment)

# generate a dummy spillover matrix
n_chs <- ncol(ss_exp)
sm <- diag(1, n_chs, n_chs)
inds <- cbind(seq_len(n_chs - 1), seq_len(n_chs)[-1])
sm[inds] <- seq(0.01, 0.3, l = n_chs - 1)
rownames(sm) <- colnames(sm) <- colnames(ss_exp)

# add artificial spillover
ss_spill <- ss_exp
exprs(ss_spill) <- exprs(ss_spill) %*% sm

# construct SCEs of reference & spilldata
ref <- t(exprs(ss_exp))
x <- prepData(ss_spill)

test_that("compCytof() - method = 'flow'/'nnls'", {
    y <- compCytof(x, sm, method = "flow", overwrite = TRUE)
    comped <- assay(y, "counts")
    expect_equivalent(ref, comped, tolerance = 1^-6)
    # test 'flow' method including comparison to 'flowCore' compensation
    comped_fC <- t(exprs(compensate(ss_spill, sm)))
    expect_equivalent(ref, comped, tolerance = 1^-6)
    expect_equivalent(ref, comped_fC, tolerance = 1^-6)
    # test 'nnls' method
    y <- compCytof(x, sm, method = "nnls", overwrite = TRUE)
    comped_nnls <- assay(y, "counts")
    expect_equivalent(ref, comped_nnls, tolerance = 1^-6)
})

test_that("compCytof() - overwrite = TRUE", {
    i <- sample(ncol(x), 1e3)
    y <- compCytof(x[, i], sm, overwrite = TRUE)
    z <- compCytof(x[, i], sm, overwrite = FALSE)
    expect_true(!"ncounts" %in% assayNames(x))
    expect_true(2*length(assays(y)) == length(assays(z)))
    expect_identical(assay(y, "counts"), assay(z, "compcounts"))
    expect_identical(assay(y, "exprs"), assay(z, "compexprs"))
})

test_that("compCytof() - cofactor = NULL", {
    cfs <- sample(10, ncol(ss_spill), TRUE)
    names(cfs) <- colnames(ss_spill)
    x <- prepData(ss_spill, cofactor = cfs)
    i <- sample(ncol(x), 1e3)
    y <- compCytof(x[, i], sm, cofactor = NULL)
    expect_identical(int_metadata(y)$cofactor, cfs)
})
