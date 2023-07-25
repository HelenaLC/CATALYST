suppressPackageStartupMessages({
    library(dplyr)
    library(SingleCellExperiment)
})
data(PBMC_fs, PBMC_panel, PBMC_md)
x0 <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
x <- cluster(x0, verbose = FALSE)
es <- assay(x, "exprs")

test_that("plotExprHeatmap() - scale = 'never'", {
    p <- plotExprHeatmap(x, scale = "never")
    replicate(10, {
        i <- sample(rownames(x), 1)
        s <- sample(levels(x$sample_id), 1)
        expect_identical(median(es[i, x$sample_id == s]), p@matrix[s, i])
    })
})

test_that("plotExprHeatmap() - scale = 'first'", {
    p <- plotExprHeatmap(x, scale = "first")
    rng <- apply(p@matrix, 2, range)
    expect_true(all(rng[1, ] >= 0))
    expect_true(all(rng[2, ] <= 1))
    expect_true(!all(rng[1, ] == 0))
    expect_true(!all(rng[2, ] == 1))
    
    p <- plotExprHeatmap(x, scale = "first", q = 0)
    y <- (es-rowMins(es))/(rowMaxs(es)-rowMins(es))
    z <- x; assay(z, "exprs") <- y
    expect_identical(t(p@matrix), .agg(z, "sample_id"))
})

test_that("plotExprHeatmap() - scale = 'last'", {
    p <- plotExprHeatmap(x, scale = "last")
    rng <- apply(p@matrix, 2, range)
    expect_true(all(rng[1, ] == 0))
    expect_true(all(rng[2, ] == 1))
    
    p <- plotExprHeatmap(x, scale = "last", q = 0)
    y <- .agg(x, "sample_id")
    y <- (y-rowMins(y))/(rowMaxs(y)-rowMins(y))
    expect_identical(t(p@matrix), y)
})

test_that("plotExprHeatmap() - row/col_clust", {
    # no clustering should retain original ordering
    p <- plotExprHeatmap(x, row_clust = FALSE, col_clust = FALSE)
    expect_identical(rownames(p@matrix), levels(x$sample_id))
    expect_identical(colnames(p@matrix), rownames(x))
})

test_that("plotExprHeatmap()", {
    # without scaling & annotations
    p <- plotExprHeatmap(x, scale = "never", bars = FALSE)
    expect_is(p, "Heatmap")
    expect_is(m <- p@matrix, "matrix")
    df <- data.frame(t(es), sample_id = x$sample_id)
    ms <- group_map(
        group_by(df, sample_id), 
        ~colMedians(as.matrix(.x), useNames=FALSE))
    ms <- do.call(rbind, ms)
    expect_identical(c(m), c(ms))
})

test_that("plotExperHeatmap() - row_anno", {
    p <- plotExprHeatmap(x, row_anno = FALSE)
    expect_is(p, "Heatmap")
    expect_true(is.null(p@left_annotation))
    
    cd_vars <- names(colData(x))
    rmv <- c("sample_id", "cluster_id")
    inc <- setdiff(cd_vars, rmv)
    
    p <- plotExprHeatmap(x, row_anno = TRUE)
    expect_is(p, "Heatmap")
    y <- p@left_annotation
    expect_is(y, "HeatmapAnnotation")
    expect_identical(names(y), inc)
    
    which <- sample(inc, 1)
    p <- plotExprHeatmap(x, row_anno = which)
    expect_is(p, "Heatmap")
    y <- p@left_annotation
    expect_is(y, "HeatmapAnnotation")
    expect_true(names(y) == which)
    
    x$foo <- sample(c("a", "b"), ncol(x), TRUE)
    p <- plotExprHeatmap(x, row_anno = TRUE)
    expect_is(p, "Heatmap")
    y <- p@left_annotation
    expect_is(y, "HeatmapAnnotation")
    expect_true(!"foo" %in% names(y))
    
    p <- plotExprHeatmap(x, row_anno = "foo")
    expect_is(p, "Heatmap")
    expect_true(is.null(p@left_annotation))
})

test_that("plotExprHeatmap() - filtering", {
    cd_vars <- names(colData(x))
    cd_vars <- setdiff(cd_vars, "cluster_id")
    # sample random factor to exclude
    # & subset of features to include
    replicate(3, {
        i <- sample(cd_vars, 1)
        ex <- sample(levels(x[[i]]), 1)
        fil <- parse(text = sprintf("%s != '%s'", i, ex))
        y <- filterSCE(x, eval(fil))
        fs <- sample(rownames(x), 5)
        p <- plotExprHeatmap(y, fs)
        expect_is(p, "Heatmap")
        expect_identical(colnames(p@matrix), fs)
        expect_identical(rownames(p@matrix), levels(y$sample_id))
    })
})

test_that("plotExprHeatmap() - by = 'cluster_id'", {
    expect_error(plotExprHeatmap(x0, by = "cluster_id"))
    k <- sample(names(cluster_codes(x))[-1], 1)
    kids <- cluster_ids(x, k)
    p <- plotExprHeatmap(x, by = "cluster_id", k = k, scale = "never")
    expect_is(p, "Heatmap")
    expect_identical(rownames(p@matrix), levels(kids))
    expect_identical(colnames(p@matrix), rownames(x))
    replicate(5, {
        k <- sample(levels(kids), 1)
        f <- sample(rownames(x), 1)
        m <- median(es[f, kids == k])
        expect_identical(p@matrix[k, f], m)
    })
})

test_that("plotExprHeatmap() - by = 'both'", {
    expect_error(plotExprHeatmap(x0, rownames(x)[1], by = "both"))
    expect_error(plotExprHeatmap(x, by = "both"))
    replicate(3, {
        f <- sample(rownames(x), 1)
        k <- sample(names(cluster_codes(x))[-1], 1)
        kids <- cluster_ids(x, k)
        p <- plotExprHeatmap(x, f, by = "both", k = k, scale = "never")
        expect_is(p, "Heatmap")
        expect_identical(rownames(p@matrix), levels(kids))
        expect_identical(colnames(p@matrix), levels(x$sample_id))
        replicate(3, {
            k <- sample(levels(kids), 1)
            s <- sample(levels(x$sample_id), 1)
            cs <- kids == k & x$sample_id == s
            m <- ifelse(!any(cs), 0, median(es[f, cs]))
            expect_identical(p@matrix[k, s], m)
        })
    })
})
