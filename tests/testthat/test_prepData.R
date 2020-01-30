context("differential")
library(flowCore)
data(PBMC_fs, PBMC_panel, PBMC_md)
chs <- PBMC_panel$fcs_colname

test_that("prepData() - cofactor", {
    expect_error(prepData(PBMC_fs, PBMC_panel, PBMC_md, cofactor = "x"))
    expect_error(prepData(PBMC_fs, PBMC_panel, PBMC_md, cofactor = c(1, 2)))
    # no transformation
    x <- prepData(PBMC_fs, PBMC_panel, PBMC_md, cofactor = NULL)
    expect_equivalent(t(assay(x)), fsApply(PBMC_fs, exprs))
    # partially unnamed cofactors should render error
    cfs <- sample(nrow(PBMC_panel))
    names(cfs) <- sample(chs); names(cfs)[8] <- NA
    expect_error(prepData(PBMC_fs, PBMC_panel, PBMC_md, cofactor = cfs))
    # unordered but named channel-specific cofactors
    cfs <- sample(nrow(PBMC_panel))
    names(cfs) <- sample(chs)
    x <- prepData(PBMC_fs, PBMC_panel, PBMC_md, cofactor = cfs)
    es <- fsApply(PBMC_fs, exprs)
    es <- vapply(chs, function(i) 
        asinh(es[, i] / cfs[i]),
        numeric(nrow(es)))
    expect_equivalent(t(assay(x)), es[, chs])
})

test_that("prepData() - panel & metadata", {
    expect_silent(prepData(PBMC_fs, as.matrix(PBMC_panel), as.matrix(PBMC_md)))
    # without marker classes
    panel <- PBMC_panel
    panel$marker_class <- NULL
    x <- prepData(PBMC_fs, panel, PBMC_md)
    expect_true(all(rowData(x)$marker_class == "none"))
})

test_that("prepData() - features", {
    # subset of features via subsetting panel & 'feature' argument
    chs_in <- sort(sample(length(chs), 5))
    x <- prepData(PBMC_fs, PBMC_panel[chs_in, ], PBMC_md)
    y <- prepData(PBMC_fs, PBMC_panel, PBMC_md, features = chs[chs_in])
    expect_identical(x, y)
})