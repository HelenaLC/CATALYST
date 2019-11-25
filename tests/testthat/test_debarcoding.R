data(sample_ff, sample_key)
ids <- rownames(sample_key)
x <- SingleCellExperiment(assays = list(counts = t(exprs(sample_ff))))
x <- assignPrelim(x, sample_key, "counts", FALSE)

test_that("assignPrelim()", {
    expect_error(assignPrelim(x, sample_key))
    expect_true(all(table(x$bc_id) == 1e3))
    expect_true(all(x$bc_id %in% rownames(sample_key)))
})

test_that("estCutoffs()", {
    x <- estCutoffs(x)
    y <- metadata(x)$sep_cutoffs
    expect_is(y, "numeric")
    expect_true(all(y[!is.na(y)] <= 1))
    expect_true(all(y[!is.na(y)] >= 0))
    expect_identical(names(y), rownames(sample_key))
    z <- x; z$bc_id <- NULL; expect_error(estCutoffs(z))
    z <- x; z$delta <- NULL; expect_error(estCutoffs(z))
})

test_that("applyCutoffs()", {
    y <- applyCutoffs(x, assay = "counts", mhl_cutoff = Inf, sep_cutoffs = 0)
    expect_identical(y$bc_id, x$bc_id)
    y <- applyCutoffs(x, assay = "counts", mhl_cutoff = 0, sep_cutoffs = 0)
    expect_true(all(y$bc_id == 0))
    y <- applyCutoffs(x, assay = "counts", sep_cutoffs = Inf)
    expect_true(all(y$bc_id == 0))
    y <- applyCutoffs(x, assay = "counts", mhl_cutoff = Inf, sep_cutoffs = 0.5)
    expect_equal(sum(y$bc_id == 0), sum(y$delta < 0.5))
    seps <- runif(nrow(sample_key)); names(seps) <- ids
    y <- applyCutoffs(x, assay = "counts", mhl_cutoff = Inf, sep_cutoffs = seps)
    ds <- split(x$delta, x$bc_id)
    expect_equal(c(table(y$bc_id)[ids]), vapply(ids, 
        function(id) sum(ds[[id]] > seps[id]), numeric(1)))
})
