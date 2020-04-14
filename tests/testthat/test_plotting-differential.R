suppressPackageStartupMessages({
    library(dplyr)
    library(SingleCellExperiment)
})
data(PBMC_fs, PBMC_panel, PBMC_md)
x0 <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
x <- cluster(x0, verbose = FALSE)
codes <- cluster_codes(x)
es <- assay(x, "exprs")

test_that("plotCounts() - prop = FALSE", {
    p <- plotCounts(x0, 
        group_by = "sample_id", 
        color_by = "condition")
    n <- p$data$value
    names(n) <- p$data$sample_id
    # output is a ggplot-object
    expect_is(p, "ggplot")
    # sample IDs are unique
    expect_true(all(table(p$data$sample_id) == 1))
    # total no. of cells & no. of cells by sample match
    expect_equal(sum(n), ncol(x0))
    expect_true(all(n == table(x0$sample_id)[names(n)]))
    # invalid arguments fail
    expect_error(plotCounts(x0, color_by = "x"))
    expect_error(plotCounts(x0, anno = "x"))
})
test_that("plotCounts() - prop = TRUE", {
    p <- plotCounts(x0, 
        prop = TRUE,
        group_by = "condition", 
        color_by = "patient_id")
    expect_is(p, "ggplot")
    # values should be proportions 
    expect_true(all(p$data$value > 0))
    expect_true(all(p$data$value < 1))
    # should sum to one for each condition
    x <- lapply(split(p$data, p$data$condition), "[[", "value")
    expect_true(all(vapply(x, sum, numeric(1)) == 1))
    ns <- table(x0$condition, x0$patient_id)
    ps <- prop.table(ns, 1)
    expect_identical(c(ps), p$data$value)
})
test_that("plotCounts() - color_by = NULL", {
    p <- plotCounts(x0, group_by = "sample_id", color_by = NULL)
    expect_is(p, "ggplot")
    # there should be no additional variables
    expect_equal(dim(p$data), c(nlevels(x0$sample_id), 2))
    # value should be number of cells per sample
    expect_identical(tabulate(x0$sample_id), p$data$value)
})

test_that("plotNRS()", {
    # include all feaures
    p <- plotNRS(x, features = NULL)
    expect_is(p, "ggplot")
    expect_is(p$data$NRS, "numeric")
    # should have one score per marker per sample
    expect_true(nrow(p$data) == nrow(x) * nlevels(x$sample_id))
    expect_true(all(table(p$data$antigen) == nlevels(x$sample_id)))
    expect_true(all(table(p$data$sample_id) == nrow(x)))
    # include type/state-markers only
    for (marker_class in c("type", "state")) {
        p <- plotNRS(x, features = marker_class)
        n <- sum(rowData(x)$marker_class == marker_class)
        expect_is(p, "ggplot")
        expect_is(p$data$NRS, "numeric")
        expect_true(nrow(p$data) == n * nlevels(x$sample_id))
        expect_true(all(table(p$data$antigen) == nlevels(x$sample_id)))
        expect_true(all(table(p$data$sample_id) == n))
    }
    # invalid arguments fail
    expect_error(plotNRS(x, features = "x"))
    expect_error(plotNRS(x, color_by = "x"))
})

test_that("pbMDS()", {
    expect_error(pbMDS(x0, color_by = "x"))
    expect_is((p <- pbMDS(x0)), "ggplot")
    expect_identical(nrow(p$data), nlevels(x0$sample_id))
    # removal of samples shouldn't cause error
    s <- sample(levels(x0$sample_id), (n <- 3))
    expect_silent(p <- pbMDS(x0[, !x$sample_id %in% s]))
    expect_equal(nrow(p$data) + 3, nlevels(x0$sample_id))
})

test_that("plotExprs()", {
    expect_error(plotExprs(x, color_by = "x"))  
    i <- sample(rownames(x), (n <- 5))
    p <- plotExprs(x, 
        features = i, 
        color_by = "condition")
    expect_is(p, "ggplot")
    expect_equal(nrow(p$data), n*ncol(x))
    expect_equivalent(levels(p$data$antigen), i)
    expect_identical(
        p$data$expression, 
        c(t(assay(x, "exprs")[i, ])))
})

test_that("plotMedExprs()", {
    # facet by antigen
    f <- sample(rownames(x), (n <- 6))
    p <- plotMedExprs(x, k, f, facet = "antigen")
    expect_is(p, "ggplot")
    expect_true(nrow(p$data) == n*nlevels(x$sample_id))
    expect_equivalent(levels(p$data$antigen), f)
    ms <- data.frame(t(es), colData(x)) %>% 
        group_by(sample_id) %>% 
        summarize_at(f, median) %>% 
        select(all_of(f)) %>% as.matrix
    expect_identical(p$data$value, c(t(ms)))
    # facet by cluster ID
    k <- sample(names(codes), 1)
    kids <- cluster_ids(x, k)
    p <- plotMedExprs(x, k, "state", facet = "cluster_id")
    expect_is(p, "ggplot")
    expect_identical(nrow(p$data), 
        length(state_markers(x))*nlevels(x$sample_id)*nlevels(kids))
    expect_equivalent(levels(p$data$antigen), state_markers(x))
})

test_that("plotExprHeatmap()", {
    # with 0-1 scaling
    expect_is(p <- plotExprHeatmap(x, scale = TRUE), "Heatmap")
    expect_is(m <- p@matrix, "matrix")
    expect_gte(min(m), 0); expect_lte(max(m), 1)
    # without scaling & all row annotations
    p <- plotExprHeatmap(x, scale = FALSE, draw_freqs = TRUE)
    expect_is(p, "Heatmap")
    expect_is(m <- p@matrix, "matrix")
    df <- data.frame(t(es), sample_id = x$sample_id)
    ms <- group_map(group_by(df, sample_id), ~colMedians(as.matrix(.x)))
    ms <- do.call(rbind, ms)
    expect_identical(c(m), c(ms))
})

test_that("plotClusterExprs()", {
    expect_error(plotClusterExprs(x0))
    expect_error(plotClusterExprs(x, k = "x"))
    kids <- cluster_ids(x, k <- sample(names(codes), 1))
    p <- plotClusterExprs(x, k, features = NULL)
    expect_is(p, "ggplot")
    expect_true(all(table(p$data$avg) == prod(dim(x))))
    expect_identical(nrow(p$data)/2, prod(dim(x)))
    expect_identical(c(t(es)), p$data$expression[p$data$avg == "no"])
})

test_that("plotAbundances()", {
    expect_error(plotAbundances(x0, k = "x"))
    expect_error(plotAbundances(x0, by = "x"))
    expect_error(plotAbundances(x, group_by = "x"))
    expect_error(plotAbundances(x, shape_by = "x"))
    kids <- cluster_ids(x, k <- sample(names(codes), 1))
    for (by in c("sample_id", "cluster_id")) {
        p <- plotAbundances(x, k, by)
        expect_is(p, "ggplot")
        expect_equal(nrow(p$data), nlevels(x$sample_id) * nlevels(kids))
        expect_true(all(table(p$data$sample_id) == nlevels(kids)))
        expect_true(all(table(p$data$cluster_id) == nlevels(x$sample_id)))
        expect_equal(p$data$value/100, c(prop.table(table(kids, x$sample_id), 2)))
    }
})

test_that("plotCodes()", {
    kids <- cluster_ids(x, k <- sample(names(codes), 1))
    expect_is(p <- plotCodes(x, k), "ggplot")
})

test_that("plotDR()", {
    dr <- reducedDim(x <- runDR(x, "UMAP", cells = (n <- 10)))
    expect_error(plotDR(x, dr = "PCA"))
    expect_error(plotDR(x, color_by = "y"))
    # colored by condition
    expect_is((p <- plotDR(x)), "ggplot")
    expect_true(!any(is.na(p$data[c("X1", "X2")])))
    expect_true(nrow(p$data) == n*nlevels(x$sample_id))
    expect_true(all(table(p$data$sample_id) == n))
    # colored by marker expression
    i <- sample(rownames(x), 1)
    p <- plotDR(x, color_by = i, scale = FALSE)
    expect_true(all(i == p$data$variable))
    expect_identical(p$data$value, es[i, !is.na(dr[, 1])])
    # colored by cluster ID
    k <- sample(names(codes), 1)
    p <- plotDR(x, color_by = k)
    expect_true(k %in% names(p$data))
    expect_identical(p$data[[k]], cluster_ids(x, k)[!is.na(dr[, 1])])
})
