suppressPackageStartupMessages({
    library(dplyr)
    library(reshape2)
    library(RColorBrewer)
})

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
df <- data.frame(t(assay(x, "exprs")), colData(x)) 
df$cluster_id <- kids <- cluster_ids(x, k)

test_that(".check_args_plotClusterHeatmap()", {
    u <- as.list(formals("plotClusterHeatmap"))
    u$x <- x
    u$k_pal <- u$m_pal <- eval(u$k_pal)
    u$hm1_pal <- u$hm2_pal <- eval(u$hm1_pal)
    expect_silent(.check_args_plotClusterHeatmap(u))
    args <- list(
        k = "x", k = c("meta2", "meta3"), m = "x", 
        hm2 = "x", hm2 = c("x", "state"),
        assay = "x", assay = rep("counts", 2))
    for (arg in c("scale", "row_anno", "row_dend", "col_dend", "draw_freqs")) {
        l <- list("x", 2, c(TRUE, FALSE))
        names(l) <- rep(arg, length(l))
        args <- c(args, l)
    }
    for (arg in c("hm1_pal", "hm2_pal", "k_pal", "m_pal")) {
        l <- list("black", c("x", "y"), sample(10))
        names(l) <- rep(arg, length(l))
        args <- c(args, l)
    }
    for (i in seq_along(args)) {
        v <- u; v[[names(args)[i]]] <- args[[i]]
        expect_error(.check_args_plotClusterHeatmap(v))
    }
})

test_that("plotClusterHeatmap() - hm2 = 'abundances'", {
    p <- plotClusterHeatmap(x, hm2 = "abundances", k = k)
    y <- p@ht_list$frequency@matrix
    expect_is(p, "HeatmapList")
    expect_true(all(rowSums(y) == 1))
    expect_equal(c(y), c(prop.table(table(kids, x$sample_id), 1)))
})

test_that("plotClusterHeatmap() - hm2 = 'state'", {
    p <- plotClusterHeatmap(x, hm2 = "state", k = k, scale = FALSE)
    y <- p@ht_list[[2]]@matrix
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
    expect_true(length(p@ht_list) == n+1)
    meds <- group_by(df, cluster_id, sample_id, 
        .drop = FALSE) %>% summarize_at(ms, median) %>% 
        replace(., is.na(.), 0) %>% {( lapply(ms, function(m) 
            acast(., cluster_id~sample_id, value.var = m)) )}
    for (i in seq_along(ms)) {
        expect_identical(p@ht_list[[i+1]]@column_title, ms[i])
        expect_identical(p@ht_list[[i+1]]@matrix, meds[[i]])
    }
})

test_that("plotClusterHeatmap() - split by patient ID", {
    p <- plotClusterHeatmap(x, k = k, split_by = "patient_id")
    expect_is(p, "list")
    expect_true(all(sapply(p, is, "Heatmap")))
    expect_identical(length(p), nlevels(x$patient_id))
})

test_that("plotClusterHeatmap() - filtering", {
    y <- mergeClusters(x, k <- "meta20", merging_table, m <- "mm")
    kids <- levels(cluster_ids(x, k))
    for (rmv in lapply(seq_len(3), sample, x = kids)) {
        fil <- filterSCE(y, !cluster_id %in% rmv, k = m)
        p <- plotClusterHeatmap(fil, k = k, m = m)
        expect_is(p, "Heatmap")
        anno <- p@left_annotation@anno_list
        expect_identical(
            rownames(p@matrix), 
            levels(cluster_ids(fil, k)))
        expect_identical(
            anno$cluster_id@color_mapping@levels,
            levels(cluster_ids(fil, k)))
        expect_equivalent(
            sort(anno$merging_id@color_mapping@levels),
            levels(cluster_ids(fil, m)))
    }
})
