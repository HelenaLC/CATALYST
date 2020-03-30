context("differential")
library(dplyr)
library(reshape2)
library(SingleCellExperiment)
data(PBMC_fs, PBMC_panel, PBMC_md)
x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)

test_that("plotClusterHeatmap() - without having run cluster()", {
    expect_error(plotClusterHeatmap(x, k = "x"))
    expect_error(plotClusterHeatmap(x, m = "x"))
})

x <- cluster(x, verbose = FALSE)
kids <- names(codes <- cluster_codes(x))
k <- sample(kids, 1)
m <- sample(setdiff(kids, k), 1)
df <- data.frame(t(assay(x)), colData(x)) 
df$cluster_id <- kids <- cluster_ids(x, k)

test_that("plotClusterHeatmap() - hm2 = 'abundances'", {
    p <- plotClusterHeatmap(x, hm2 = "abundances", k = k)
    y <- p@ht_list$frequency@matrix
    expect_is(p, "HeatmapList")
    expect_true(all(rowSums(y) == 1))
    expect_equal(c(y), c(prop.table(table(kids, x$sample_id), 1)))
})

test_that("plotClusterHeatmap() - hm2 = 'state_markers'", {
    p <- plotClusterHeatmap(x, hm2 = "state_markers", k = k, scale = FALSE)
    y <- p@ht_list[[3]]@matrix
    expect_is(p, "HeatmapList")
    # check dimension names
    expect_identical(dim(y), c(nlevels(kids), length(state_markers(x))))
    expect_identical(rownames(y), levels(kids))
    expect_identical(colnames(y), state_markers(x))
    # check 2nd heatmap data
    ms <- group_by(df, cluster_id) %>% 
        summarize_at(state_markers(x), median) %>% 
        do(.[, state_markers(x)]) %>% do.call(what = cbind)
    expect_identical(c(y), c(ms))
})

test_that("plotClusterHeatmap() - hm2 = specific state markers", {
    ms <- sample(state_markers(x), (n <- 3))
    expect_error(plotClusterHeatmap(x, hm2 = c(ms, "x"), k = k))
    p <- plotClusterHeatmap(x, hm2 = ms, k = k, scale = FALSE)
    expect_is(p, "HeatmapList")
    expect_true(length(p@ht_list) == n+2)
    meds <- group_by(df, cluster_id, sample_id, 
        .drop = FALSE) %>% summarize_at(ms, median) %>% 
        replace(., is.na(.), 0) %>% {( lapply(ms, function(m) 
            acast(., cluster_id~sample_id, value.var = m)) )}
    for (i in seq_along(ms)) {
        expect_identical(p@ht_list[[i+2]]@column_title, ms[i])
        expect_identical(p@ht_list[[i+2]]@matrix, meds[[i]])
    }
})

test_that("plotClusterHeatmap() - split by patient ID", {
    p <- plotClusterHeatmap(x, k = k, split_by = "patient_id")
    expect_is(p, "list")
    expect_true(all(sapply(p, is, "HeatmapList")))
    expect_identical(length(p), nlevels(x$patient_id))
})

test_that("plotClusterHeatmap() - with all row annotations", {
    p <- plotClusterHeatmap(x, k = k, m = m, draw_freqs = TRUE, scale = FALSE)
    expect_is(p, "HeatmapList")
    y <- p@ht_list$expression@matrix
    # check dimension names
    expect_identical(dim(y), c(nlevels(kids), length(type_markers(x))))
    expect_identical(rownames(y), levels(kids))
    expect_identical(colnames(y), type_markers(x))
    # check 1st heatmap data
    ms <- group_by(df, cluster_id) %>% 
        summarize_at(type_markers(x), median) %>% 
        do(.[, type_markers(x)]) %>% do.call(what = cbind)
    expect_identical(c(y), c(ms))
    # check left row annotations of cluster & merging IDs
    expect_identical(c(p@ht_list$cluster_id@matrix), levels(kids))
    expect_identical(c(p@ht_list$merging_id@matrix), 
        as.character(codes[match(levels(kids), codes[[k]]), m]))
    # check right row annotation of frequency histograms
    u <- p@ht_list[[4]]@anno_list[[1]]@fun@data_scale
    expect_gte(min(u), 0); expect_equal(max(u),
        max(table(kids)/ncol(x))*105, tol = 1e-3)
    expect_equal(p@ht_list[[5]]@anno_list[[1]]@fun@data_scale, c(0, 1))
})
