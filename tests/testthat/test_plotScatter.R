data(sample_ff, sample_key)
x <- prepData(sample_ff)
# sample only a couple IDs for testing
ids <- sample(rownames(sample_key), 3)
x <- assignPrelim(x, sample_key[ids, ], verbose = FALSE)
x <- x[, x$bc_id %in% ids]

test_that("plotScatter() - labels", {
    chs <- channels(x)
    expect_error(plotScatter(x, sample(chs, 2), label = "x"))
    args <- eval(formals("plotScatter")$label)
    for (l in args) {
        ls <- switch(l, 
            target = rownames(x), channel = chs, 
            paste(chs, rownames(x), sep = "-"))
        i <- sample(nrow(x), 2)
        p <- plotScatter(x, rownames(x)[i], label = l)
        expect_is(p, "ggplot")
        labs <- unlist(p$labels[c("x", "y")])
        expect_equivalent(labs, ls[i])
    }
})

test_that("plotScatter() - facetting", {
    chs <- sample(rownames(x), 3)
    p <- plotScatter(x, chs, label = "target")
    expect_is(p, "ggplot")
    expect_true(p$labels$x == chs[1])
    expect_identical(levels(p$data$variable), chs[-1])
})

test_that("plotScatter() - color_by", {
    chs <- sample(rownames(x), 2)
    p <- plotScatter(x, chs, color_by = "delta")
    expect_is(p, "ggplot")
    expect_is(p$scales$scales[[1]], "ScaleContinuous")    
    p <- plotScatter(x, chs, color_by = "bc_id")
    expect_is(p, "ggplot")
    expect_is(p$guides$guides[[1]], "GuideLegend")
})
