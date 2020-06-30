suppressPackageStartupMessages({
    library(flowCore)
})

data(PBMC_fs, PBMC_panel, PBMC_md)
chs <- PBMC_panel$fcs_colname

test_that("prepData() - cofactor", {
    expect_error(prepData(PBMC_fs, PBMC_panel, PBMC_md, cofactor = "x"))
    expect_error(prepData(PBMC_fs, PBMC_panel, PBMC_md, cofactor = c(1, 2)))
    # no transformation
    x <- prepData(PBMC_fs, PBMC_panel, PBMC_md, transform = FALSE)
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
    expect_equivalent(t(assay(x, "exprs")), es[, chs])
})

test_that("prepData() - panel & metadata", {
    expect_silent(prepData(PBMC_fs, as.matrix(PBMC_panel), as.matrix(PBMC_md)))
    # with invalid marker classes
    panel <- PBMC_panel
    panel$marker_class <- "x"
    expect_error(prepData(PBMC_fs, panel, PBMC_md))
    # without marker classes
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

test_that("prepData() - panel & md = NULL", {
    # construct artificial flowFrame
    n <- 100; m <- 50
    x <- matrix(runif(n * m), n, m)
    colnames(x) <- seq_len(m)
    ff <- flowFrame(x)
    # dimensions should be reversed
    y <- prepData(ff, by_time = FALSE)
    expect_is(y, "SingleCellExperiment")
    expect_identical(rev(dim(y)), dim(x))
    # data should be unchanged
    expect_equal(t(assay(y)), x, tolerance = 1e-6)
    # with transformation
    expect_error(prepData(ff, transform = "x"))
    y <- prepData(ff, by_time = FALSE,
        transform = TRUE, cofactor = (cf <- 20))
    expect_true("exprs" %in% assayNames(y))
    expect_equivalent(assay(y, "exprs"), asinh(t(exprs(ff))/cf))
    # construct artifical flowSet
    i <- sample(seq_len(n), 10)
    j <- seq_len(n)[-i]
    fs <- flowSet(ff[i, ], ff[j, ])
    # should render message when 
    # acquisition times aren't available
    expect_message(prepData(fs))
    # assure flowFrames are ordered correctly
    keyword(fs[[1]])$`$BTIM` <- 2
    keyword(fs[[2]])$`$BTIM` <- 1
    expect_silent(y <- prepData(fs))
    expect_equal(assay(y)[, seq_along(j)], t(exprs(ff[j, ])))
    expect_equal(assay(y)[, seq_along(i)+length(j)], t(exprs(ff[i, ])))
})

test_that("prepData() - panel_cols", {
    # all panel_cols differ from default
    panel <- PBMC_panel
    panel_cols <- list(channel = "channel", antigen = "target", class = "class")
    names(panel) <- unlist(panel_cols)
    x <- prepData(PBMC_fs, panel, PBMC_md, panel_cols = panel_cols)
    y <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
    expect_identical(x, y)
    # some panel_cols differ from default
    panel <- PBMC_panel
    panel_cols <- list(channel = "channel")
    names(panel)[1] <- unlist(panel_cols)
    x <- prepData(PBMC_fs, panel, PBMC_md, panel_cols = panel_cols)
    y <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
    expect_identical(x, y)
})

test_that("prepData() - md_cols", {
    # all md_cols differ from default
    md_cols <- list(file = "file", id = "id", factors = c("group", "patient"))
    md <- PBMC_md; names(md) <- unlist(md_cols)
    x <- prepData(PBMC_fs, PBMC_panel, md, md_cols = md_cols)
    y <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
    expect_equivalent(x, y)
    # some md_cols differ from default
    md_cols <- list(file = "file", id = "id")
    md <- PBMC_md; names(md)[c(1, 2)] <- unlist(md_cols)
    x <- prepData(PBMC_fs, PBMC_panel, md, md_cols = md_cols)
    y <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
    expect_identical(x, y)
})
