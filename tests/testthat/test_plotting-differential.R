context("differential")
library(dplyr)
library(SingleCellExperiment)
data(PBMC_fs, PBMC_panel, PBMC_md)
x0 <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
x <- cluster(x0, verbose = FALSE)
codes <- cluster_codes(x)

test_that("plotCounts()", {
    p <- plotCounts(x0, color_by = "condition")
    n <- p$data$n_cells
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

test_that("plotMDS()", {
    expect_error(plotMDS(x0, color_by = "x"))
    expect_is((p <- plotMDS(x0)), "ggplot")
    expect_identical(nrow(p$data), nlevels(x0$sample_id))
    # removal of samples shouldn't cause error
    s <- sample(levels(x0$sample_id), (n <- 3))
    expect_silent(p <- plotMDS(x0[, !x$sample_id %in% s]))
    expect_equal(nrow(p$data) + 3, nlevels(x0$sample_id))
})

test_that("plotExprs()", {
    expect_error(plotExprs(x0, color_by = "x"))  
    p <- plotExprs(x0, color_by = "condition")
    expect_is(p, "ggplot")
    expect_equal(nrow(p$data), prod(dim(x0)))
    expect_identical(p$data$expression, c(t(assay(x0))))
})

test_that("plotMedExprs()", {
    k <- sample(names(codes), 1)
    kids <- cluster_ids(x, k)
    
    # facet by antigen
    p <- plotMedExprs(x, k, facet = "antigen")
    expect_is(p, "ggplot")
    expect_identical(nrow(p$data), nrow(x)*nlevels(x$sample_id))
    ms <- data.frame(t(assay(x)), colData(x)) %>% 
        group_by(sample_id) %>% 
        summarize_at(rownames(x), median) %>% 
        select(rownames(x)) %>% as.matrix
    expect_identical(p$data$value, c(t(ms)))
    
    # facet by cluster ID
    p <- plotMedExprs(x, k, facet = "cluster_id")
    n <- sum(rowData(x)$marker_class == "state")
    expect_is(p, "ggplot")
    expect_identical(nrow(p$data), n*nlevels(x$sample_id)*nlevels(kids))
})

test_that("plotExprHeatmap()", {
    # with 0-1 scaling
    expect_is(p <- plotExprHeatmap(x, scale = TRUE), "HeatmapList")
    expect_is(m <- p@ht_list$expression@matrix, "matrix")
    expect_gte(min(m), 0); expect_lte(max(m), 1)
    # without scaling & all row annotations
    expect_is(p <- plotExprHeatmap(x, scale = FALSE, draw_freqs = TRUE), "HeatmapList")
    expect_is(m <- p@ht_list$expression@matrix, "matrix")
    expect_is(anno <- p@ht_list[[3]], "HeatmapAnnotation")
    rng <- anno@anno_list[[1]]@fun@data_scale
    expect_equal(rng, c(0, max(table(x$sample_id)*1.05)))
    df <- data.frame(t(assay(x)), id = colData(x)$sample_id)
    ms <- group_map(group_by(df, id), ~colMedians(as.matrix(.x)))
    ms <- do.call(rbind, ms)
    expect_identical(c(m), c(ms))
})

test_that("plotClusterExprs()", {
    expect_error(plotClusterExprs(x0))
    expect_error(plotClusterExprs(x, k = "x"))
    kids <- cluster_ids(x, k <- sample(names(codes), 1))
    expect_is(p <- plotClusterExprs(x, k), "ggplot")
    expect_true(all(table(p$data$ref) == prod(dim(x))))
    expect_identical(nrow(p$data)/2, prod(dim(x)))
    expect_identical(c(t(assay(x))), p$data$expression[p$data$ref == "no"])
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
        expect_equal(p$data$freq/100, c(prop.table(table(kids, x$sample_id), 2)))
    }
})

test_that("plotCodes()", {
    kids <- cluster_ids(x, k <- sample(names(codes), 1))
    expect_is(p <- plotCodes(x, k, verbose = FALSE), "ggplot")
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
    p <- plotDR(x, color_by = i)
    expect_true(i %in% names(p$data))
    expect_identical(p$data[[i]], assay(x)[i, !is.na(dr[, 1])])
    # colored by cluster ID
    k <- sample(names(codes), 1)
    p <- plotDR(x, color_by = k)
    expect_true(k %in% names(p$data))
    expect_identical(p$data[[k]], cluster_ids(x, k)[!is.na(dr[, 1])])
})
