suppressPackageStartupMessages({
    library(flowCore)
    library(SingleCellExperiment)
})

x <- .toySCE()
es <- assay(x, "exprs")

test_that(".check_colors()", {  
    expect_error(.check_colors("blue")) 
    expect_error(.check_colors(letters))
    expect_silent(.check_colors(colors()))
    expect_error(.check_colors(c("x", colors())))
})

test_that("guessPanel()", {
    expect_error(guessPanel("x"))
    ff <- get(data("sample_ff"))
    ps0 <- parameters(ff)
    # channel & antigen in 'name' & 'desc'
    expect_is(df <- guessPanel(ff), "data.frame")
    expect_true(nrow(df) == ncol(ff))
    expect_true(all(df$fcs_colname == ps0@data$name))
    expect_true(all(df$antigen == ps0@data$desc))
    expect_true(all(df$desc == ps0@data$desc))
    # 'desc' of the from channel_antigen
    ps <- ps0
    ps@data$desc <- with(ps@data, paste(name, desc, sep = "_"))
    parameters(ff) <- ps
    expect_is(df <- guessPanel(ff), "data.frame")
    expect_true(all(df$fcs_colname == ps@data$name))
    expect_true(all(df$antigen == ps0@data$desc))
    expect_true(all(df$desc == ps0@data$name))
})

test_that(".split_cells() by 1 factor", {
    for (by in c("sample_id", "cluster_id")) {
        cs <- .split_cells(x, by)
        expect_identical(names(cs), levels(x[[by]]))
        expect_identical(
            as.numeric(vapply(cs, length, numeric(1))),
            as.numeric(table(x[[by]])))
    }
})

test_that(".split_cells() by 2 factors", {
    for (by in list(
        c("sample_id", "cluster_id"), 
        c("cluster_id", "sample_id"))) {
        cs <- .split_cells(x, by)
        expect_identical(names(cs), levels(x[[by[1]]]))
        expect_true(all(vapply(cs, function(u) 
            identical(names(u), levels(x[[by[2]]])), 
            logical(1))))
    }
})

test_that(".agg() by 1 factor", {
    for (by in c("sample_id", "cluster_id")) {
        for (fun in c("sum", "mean", "median")) {
            pb <- .agg(x, by, fun)
            replicate(10, {
                i <- sample(seq_len(nrow(x)), 1)
                j <- sample(levels(x[[by]]), 1)
                expect_equal(pb[i, j], get(fun)(es[i, x[[by]] == j]))
            })
        }
    }
})

test_that(".agg() by 2 factors", {
    for (by in list(
        c("sample_id", "cluster_id"), 
        c("cluster_id", "sample_id"))) {
        for (fun in c("sum", "mean", "median")) {
            pb <- .agg(x, by, fun, "exprs")
            expect_identical(names(pb), levels(x[[by[1]]]))
            expect_true(all(vapply(pb, function(u) 
                identical(colnames(u), levels(x[[by[2]]])), 
                logical(1))))
            replicate(10, {
                a <- sample(levels(x[[by[1]]]), 1)
                b <- sample(levels(x[[by[2]]]), 1)
                i <- sample(nrow(x), 1)
                j <- x[[by[1]]] == a & x[[by[2]]] == b
                expect_equal(pb[[a]][i, b], get(fun)(es[i, j]))
            })
        }
    }
})
