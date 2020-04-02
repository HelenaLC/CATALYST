context("utilities for differential analysis")
set.seed(as.numeric(format(Sys.time(), "%s")))
x <- .toySCE()
y <- assay(x, "exprs")

test_that(".get_features", {
    mcs <- sample(c("type", "state"), nrow(x), TRUE)
    rowData(x)$marker_class <- mcs
    expect_identical(.get_features(x, NULL), rownames(x))
    for (class in c("type", "state"))
        expect_identical(
            .get_features(x, class), 
            get(paste0(class, "_markers"))(x))
    fs <- sample(rownames(x), 10)
    expect_identical(.get_features(x, fs), fs)
    fs[sample(10, 1)] <- "x"
    expect_error(.get_features(x, "x"))
    expect_error(.get_features(x, fs))
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
                expect_equal(pb[i, j], get(fun)(y[i, x[[by]] == j]))
            })
        }
    }
})

test_that(".agg() by 2 factors", {
    for (by in list(
        c("sample_id", "cluster_id"), 
        c("cluster_id", "sample_id"))) {
        for (fun in c("sum", "mean", "median")) {
            pb <- .agg(x, by, fun)
            expect_identical(names(pb), levels(x[[by[1]]]))
            expect_true(all(vapply(pb, function(u) 
                identical(colnames(u), levels(x[[by[2]]])), 
                logical(1))))
            replicate(10, {
                a <- sample(levels(x[[by[1]]]), 1)
                b <- sample(levels(x[[by[2]]]), 1)
                i <- sample(seq_len(nrow(x)), 1)
                j <- x[[by[1]]] == a & x[[by[2]]] == b
                expect_equal(pb[[a]][i, b], get(fun)(y[i, j]))
            })
        }
    }
})
