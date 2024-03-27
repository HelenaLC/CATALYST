suppressPackageStartupMessages({
    library(dplyr)
    library(reshape2)
    library(RColorBrewer)
})

data(PBMC_fs, PBMC_panel, PBMC_md)
x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)

test_that("plotMultiHeatmap() - without having run cluster()", {
    expect_error(plotMultiHeatmap(x, k = "x"))
    expect_error(plotMultiHeatmap(x, m = "x"))
})

x <- cluster(x, verbose = FALSE)
kids <- names(codes <- cluster_codes(x))
k <- sample(kids, 1)
m <- sample(setdiff(kids, k), 1)
df <- data.frame(t(assay(x, "exprs")), colData(x)) 
df$cluster_id <- kids <- cluster_ids(x, k)

test_that("plotMultiHeatmap() - hm2 = 'abundances'", {
    p <- plotMultiHeatmap(x, hm2 = "abundances", k = k, normalize = FALSE)
    y <- p@ht_list$frequency@matrix
    expect_is(p, "HeatmapList")
    expect_equivalent(colSums(y), rep(1, ncol(y)))
    expect_equal(c(y), c(prop.table(table(kids, x$sample_id), 2)))
})

test_that("plotMultiHeatmap() - hm2 = 'state'", {
    p <- plotMultiHeatmap(x, hm2 = "state", k = k, scale = "never")
    y <- p@ht_list[[2]]@matrix
    expect_is(p, "HeatmapList")
    # check dimension names
    expect_identical(dim(y), c(nlevels(kids), length(state_markers(x))))
    expect_identical(rownames(y), levels(kids))
    expect_identical(colnames(y), state_markers(x))
    # check 2nd heatmap data
    ms <- group_by(df, cluster_id) %>% 
        summarize_at(state_markers(x), median) %>% 
        do(.[, state_markers(x)]) %>% 
        do.call(what = "cbind")
    expect_equal(c(y), c(ms))
})

test_that("plotMultiHeatmap() - hm2 = specific state markers", {
    ms <- sample(state_markers(x), (n <- 3))
    expect_error(plotMultiHeatmap(x, hm2 = c(ms, "x"), k = k))
    p <- plotMultiHeatmap(x, hm2 = ms, k = k, scale = "never")
    expect_is(p, "HeatmapList")
    expect_true(length(p@ht_list) == n+1)
    meds <- group_by(df, cluster_id, sample_id, 
        .drop = FALSE) %>% summarize_at(ms, median) %>% 
        replace(., is.na(.), 0) %>% {( lapply(ms, function(m) 
            acast(., cluster_id~sample_id, value.var = m)) )}
    for (i in seq_along(ms)) {
        expect_equal(p@ht_list[[i+1]]@column_title, ms[i])
        expect_equal(p@ht_list[[i+1]]@matrix, meds[[i]])
    }
})

test_that("plotMultiHeatmap() - filtering", {
    y <- mergeClusters(x, k <- "meta20", merging_table, m <- "mm")
    kids <- levels(cluster_ids(x, k))
    for (rmv in lapply(seq_len(3), sample, x = kids)) {
        fil <- filterSCE(y, !cluster_id %in% rmv, k = m)
        p <- plotMultiHeatmap(fil, k = k, m = m)
        expect_is(p, "HeatmapList")
        anno <- p@ht_list[[1]]@left_annotation@anno_list
        expect_identical(
            rownames(p@ht_list[[1]]@matrix), 
            levels(cluster_ids(fil, k)))
        expect_identical(
            anno$cluster_id@color_mapping@levels,
            levels(cluster_ids(fil, k)))
        expect_equivalent(
            sort(anno$merging_id@color_mapping@levels),
            levels(cluster_ids(fil, m)))
    }
})
