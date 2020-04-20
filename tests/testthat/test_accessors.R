suppressPackageStartupMessages({
    library(SingleCellExperiment)
})
data(PBMC_fs, PBMC_panel, PBMC_md)
x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
x <- cluster(x, verbose = FALSE)

test_that("ei()", {
    y <- x; class(y) <- "x"
    expect_error(ei(y))
    expect_identical(ei(x), metadata(x)$experiment_info)
    u <- ei(x)
    for (i in colnames(u)[-ncol(u)]) {
        j <- sample(levels(u[[i]]), 1)
        v <- u[u[[i]] != j, ]
        fil <- parse(text = sprintf("%s != '%s'", i, j))
        expect_equivalent(ei(filterSCE(x, eval(fil))), v)
    }
})

test_that("n_cells()", {
    y <- x; class(y) <- "x"
    expect_error(n_cells(y))
    ns <- n_cells(x)
    expect_identical(sum(ns), ncol(x))
    expect_identical(c(unname(ns)), tabulate(x$sample_id))
    sids <- sample(levels(x$sample_id), 3)
    cs <- x$sample_id %in% sids
    x <- x[, cs]
    ns <- n_cells(x)
    expect_true(sum(ns) == sum(cs))
    expect_true(setequal(names(ns), sids))
    x$sample_id <- NULL
    expect_silent(ns <- n_cells(x))
    expect_true(is.null(ns))
})

test_that("channels()", {
    y <- x; class(y) <- "x"
    expect_error(channels(y))
    expect_equivalent(channels(x), rowData(x)$channel_name)
    rowData(x)$channel_name <- NULL
    expect_silent(channels(x))
})

test_that("marker_classes()", {
    y <- x; class(y) <- "x"
    expect_error(marker_classes(y))
    expect_equivalent(marker_classes(x), rowData(x)$marker_class)
    rowData(x)$marker_class <- NULL
    expect_silent(marker_classes(x))
})

test_that("type_markers()", {
    y <- x; class(y) <- "x"
    expect_error(type_markers(y))
    expect_equivalent(type_markers(x),
        rownames(x)[marker_classes(x) == "type"])
    y <- x; rowData(y)$marker_class <- NULL
    expect_silent(type_markers(y))
    rowData(x)$marker_class <- "x"
    expect_silent(type_markers(x))
})

test_that("state_markers()", {
    y <- x; class(y) <- "x"
    expect_error(state_markers(y))
    expect_equivalent(state_markers(x),
        rownames(x)[marker_classes(x) == "state"])
    y <- x; rowData(y)$marker_class <- NULL
    expect_silent(state_markers(y))
    rowData(x)$marker_class <- "x"
    expect_silent(state_markers(x))
})

test_that("cluster_ids()", {
    y <- x; class(y) <- "x"
    expect_error(cluster_ids(y))
    kids <- names(cluster_codes(x))
    expect_identical(cluster_ids(x), x$cluster_id)
    expect_identical(cluster_ids(x, kids[1]), x$cluster_id)
    expect_error(cluster_ids(x, "x"))
    ks <- sample(kids[-1], 5)
    expect_error(cluster_ids(x, ks))
    kids <- as.numeric(as.character(x$cluster_id))
    for (k in ks)
        expect_identical(
            cluster_ids(x, k), 
            cluster_codes(x)[[k]][kids])
})

test_that("delta_area()", {
    y <- x; class(y) <- "x"
    expect_error(delta_area(y))
    expect_is(p <- delta_area(x), "ggplot")
    expect_identical(p, metadata(x)$delta_area)
    metadata(x)$delta_area <- NULL
    expect_silent(delta_area(x))
})
