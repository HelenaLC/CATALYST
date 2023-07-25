suppressPackageStartupMessages({
    library(dplyr)
    library(SingleCellExperiment)
})
set.seed(3004)
data(PBMC_fs, PBMC_panel, PBMC_md)
x0 <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
x <- cluster(x0, seed = 3004, verbose = FALSE)
es <- assay(x, "exprs")
codes <- cluster_codes(x)

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

test_that("pbMDS() - by = 'sample_id'", {
    expect_error(pbMDS(x0, color_by = "x"))
    expect_is((p <- pbMDS(x0)), "ggplot")
    expect_identical(nrow(p$data), nlevels(x0$sample_id))
    # removal of samples
    s <- sample(levels(x0$sample_id), (n <- 3))
    expect_is(p <- pbMDS(filterSCE(x0, !sample_id %in% s)), "ggplot")
    expect_equal(nrow(p$data) + 3, nlevels(x0$sample_id))
})
test_that("pbMDS() - by = 'cluster_id'", {
    k <- sample(names(codes)[-seq_len(5)], 1)
    nk <- length(kids <- levels(codes[[k]]))
    expect_is(p <- pbMDS(x, by = "cluster_id", k = k), "ggplot")
    expect_identical(nrow(p$data), nk)
    expect_identical(levels(p$data$cluster_id), kids)
    expect_equivalent(p$data$n_cells, c(table(cluster_ids(x, k))))
    # removal of clusters
    ks <- sample(kids, 3)
    y <- filterSCE(x, !cluster_id %in% ks, k = k)
    expect_is(p <- pbMDS(y, by = "cluster_id", k = k), "ggplot")
    expect_identical(levels(p$data$cluster_id), setdiff(kids, ks))
    expect_equivalent(p$data$n_cells, c(table(cluster_ids(y, k))))
})
test_that("pbMDS() - by = 'both'", {
    k <- sample(names(codes)[-seq_len(5)], 1)
    nk <- length(kids <- levels(codes[[k]]))
    ns <- length(sids <- levels(x$sample_id))
    ks <- sample(kids, 3); ss <- sample(sids, 3)
    .check_nc <- function() {
        nc <- table(cluster_ids(y, k), y$sample_id)
        nc <- nc[as.matrix(p$data[c("cluster_id", "sample_id")])]
        expect_identical(p$data$n_cells, nc)
    }
    y <- filterSCE(x, k = k, !(cluster_id %in% ks & sample_id %in% ss)) 
    expect_is(p <- pbMDS(y, by = "both", k = k), "ggplot")
    expect_identical(nrow(p$data), nk*ns)
    expect_true(all(table(p$data$sample_id) == nk))
    expect_true(all(table(p$data$cluster_id) == ns))
    .check_nc()
    # removal of clusters 
    y <- filterSCE(x, k = k, !cluster_id %in% ks) 
    expect_is(p <- pbMDS(y, by = "both", k = k), "ggplot")
    expect_identical(levels(p$data$cluster_id), setdiff(kids, ks))
    .check_nc()
    nc <- table(cluster_ids(y, k), y$sample_id)
    nc <- nc[as.matrix(p$data[c("cluster_id", "sample_id")])]
    expect_identical(p$data$n_cells, nc)
    # removal of samples 
    y <- filterSCE(x, k = k, !sample_id %in% ss) 
    expect_is(p <- pbMDS(y, by = "both", k = k), "ggplot")
    expect_identical(levels(p$data$sample_id), setdiff(sids, ss))
    .check_nc()
})
test_that("pbMDS() - single factor in colData", {
    colData(x) <- colData(x)[, seq(2)]
    p <- pbMDS(x, by = "sample_id")
    p <- tryCatch(print(p), error = function(e) e)
    expect_false(inherits(p, "error"))
})

test_that("clrDR()", {
    k <- sample(names(codes)[-seq_len(11)], 1)
    nk <- nlevels(kids <- cluster_ids(x, k))
    ns <- nlevels(sids <- x$sample_id)
    # by = 'sample_id'
    expect_is(p <- clrDR(x, by = "sample_id", k = k), "ggplot")
    expect_identical(nrow(p$data), nlevels(x$sample_id))
    expect_identical(p$data$n_cells, tabulate(x$sample_id))
    # PC loading arrows
    ls <- p$layers; y <- ls[[length(ls)]]
    expect_identical(nrow(y$data), nlevels(kids))
    expect_equivalent(levels(y$data$cluster_id), levels(kids))
    # removal of samples
    ss <- sample(levels(sids), n <- 3)
    y <- filterSCE(x, k = k, !sample_id %in% ss) 
    p <- clrDR(y, by = "sample_id", k = k)
    expect_identical(nrow(p$data), as.integer(ns-n))
    expect_identical(p$data$n_cells, tabulate(y$sample_id))
    expect_identical(levels(p$data$sample_id), setdiff(levels(sids), ss))
    # by = 'cluster_id'
    expect_is(p <- clrDR(x, by = "cluster_id", k = k), "ggplot")
    expect_identical(nrow(p$data), nk)
    expect_identical(p$data$n_cells, tabulate(kids))
    expect_identical(levels(p$data$cluster_id), levels(kids))
    # PC loading arrows
    ls <- p$layers; y <- ls[[length(ls)]]
    expect_identical(nrow(y$data), ns)
    expect_equivalent(levels(y$data$sample_id), levels(sids))
    # removal of clusters
    ks <- sample(levels(kids), n <- 3)
    y <- filterSCE(x, k = k, !cluster_id %in% ks) 
    p <- clrDR(y, by = "cluster_id", k = k)
    expect_identical(nrow(p$data), as.integer(nk-n))
    expect_identical(p$data$n_cells, tabulate(cluster_ids(y, k)))
    expect_identical(levels(p$data$cluster_id), setdiff(levels(kids), ks))
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

test_that("plotPbExprs()", {
    # facet by antigen
    f <- sample(rownames(x), (n <- 6))
    p <- plotPbExprs(x, k, f, facet = "antigen")
    expect_is(p, "ggplot")
    expect_true(nrow(p$data) == n*nlevels(x$sample_id))
    expect_equivalent(levels(p$data$antigen), f)
    ms <- data.frame(t(es), colData(x), check.names = FALSE) %>% 
        group_by(sample_id) %>% 
        summarize_at(f, median) %>% 
        select(all_of(f)) %>% as.matrix
    expect_identical(p$data$value, c(t(ms)))
    # facet by cluster ID
    k <- sample(names(codes), 1)
    kids <- cluster_ids(x, k)
    p <- plotPbExprs(x, k, "state", facet = "cluster_id")
    expect_is(p, "ggplot")
    expect_identical(nrow(p$data), 
        length(state_markers(x))*sum(table(x$sample_id, kids) != 0))
    expect_equivalent(levels(p$data$antigen), state_markers(x))
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
        expect_equal(p$data$Freq/100, c(prop.table(table(kids, x$sample_id), 2)))
    }
})

test_that("plotAbundances() - filtering", {
    # exclude random subset of clusters
    k <- sample(names(codes)[-seq_len(4)], 1)
    kids <- levels(cluster_ids(x, k))
    ns <- sample(length(kids)-1, 3)
    for (ks in lapply(ns, sample, x = kids)) {
        y <- filterSCE(x, !cluster_id %in% ks, k = k)
        p <- plotAbundances(y, k, by = "sample_id")
        expect_is(p, "ggplot")
        expect_true(setequal(p$data$cluster_id, setdiff(kids, ks)))
    }
    # exclude random subset of samples
    sids <- levels(x$sample_id)
    ns <- sample(length(sids)-1, 3)
    for (ss in lapply(ns, sample, x = sids)) {
        y <- filterSCE(x, !sample_id %in% ss, k = k)
        p <- plotAbundances(y, k, by = "sample_id")
        expect_is(p, "ggplot")
        expect_true(setequal(levels(p$data$sample_id), setdiff(sids, ss)))
    }
})

test_that("plotCodes()", {
    kids <- cluster_ids(x, k <- sample(names(codes), 1))
    expect_is(p <- plotCodes(x, k), "ggplot")
})

test_that("plotDR()", {
    dr <- reducedDim(x <- runDR(x, "UMAP", cells = (n <- 10)))
    expect_error(plotDR(x, dr = "x"))
    expect_error(plotDR(x, color_by = "x"))
    # colored by condition
    expect_is(p <- plotDR(x), "ggplot")
    expect_true(!any(is.na(p$data[c("x", "y")])))
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
