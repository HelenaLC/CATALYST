data(PBMC_fs, PBMC_panel, PBMC_md)
x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
x <- x[, sample(ncol(x), 500)]

test_that("cluster() - default parameters", {
    expect_error(cluster(x, features = "x"))
    set.seed(1); x1 <- cluster(x, maxK = 3, verbose = FALSE)
    set.seed(1); x2 <- cluster(x, maxK = 3, 
        features = type_markers(x), verbose = FALSE)
    expect_identical(x1$cluster_id, x2$cluster_id)
})

test_that("cluster() - use feature-subset", {
    y <- x[(fs <- sample(rownames(x), 5)), ]
    rowData(y)$marker_class <- "type"
    set.seed(1); x1 <- cluster(y, maxK = 3, verbose = FALSE)
    set.seed(1); x2 <- cluster(x, maxK = 3, features = fs, verbose = FALSE)
    expect_identical(x1$cluster_id, x2$cluster_id)
})

kid <- paste0("meta", (k <- 10))
x <- cluster(x, maxK = k, verbose = FALSE)
tbl <- data.frame(seq_len(k), sample(seq_len(k), k, replace = TRUE))
tbl <- dplyr::mutate_all(tbl, as.character)
    
test_that("mergeClusters()", {
    y <- mergeClusters(x, kid, tbl, (id <- "id"))
    expect_error(mergeClusters(y, kid, tbl, "id"))
    expect_identical(names(cluster_codes(y))[k+1], id)
    expect_identical(cluster_codes(y)[-(k+1)], cluster_codes(x))
    ns <- sapply(sort(unique(tbl[[2]])), function(id)
        sum(cluster_ids(x, kid) %in% tbl[[1]][tbl[[2]] == id]))
    expect_identical(c(table(cluster_ids(y, "id"))), ns)
})
