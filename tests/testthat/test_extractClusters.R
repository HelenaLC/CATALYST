sce0 <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
sce <- cluster(sce0, verbose = FALSE)

test_that("extractClusters() - without having run cluster()", {
    expect_error(extractClusters(sce0))
})

test_that("extractClusters()", {
    # check vebosity
    expect_message(extractClusters(sce, NULL))
    expect_silent(extractClusters(sce, NULL, verbose = FALSE))
    # check fails when specifying invalid clustering or clusters
    expect_error(extractClusters(sce, "x", verbose = FALSE))
    expect_error(extractClusters(sce, "meta10", 11, verbose = FALSE))
    # check default
    x <- extractClusters(sce, NULL, verbose = FALSE)
    expect_is(x, "flowSet")
    expect_true(all(unlist(fsApply(x, function(u) is(u, "flowFrame")))))
    expect_identical(c(fsApply(x, identifier)), levels(cluster_ids(sce)))
    expect_identical(
        as.numeric(fsApply(x, nrow)), 
        as.numeric(table(cluster_ids(sce))))
    # 5 random spot-checks
    replicate(5, {
        # sample clustering & clusters to extract
        k <- sample(names(cluster_codes(sce)), 1)
        kids <- cluster_ids(sce, k)
        nk <- sample(nlevels(kids), 1)
        ks <- sample(levels(kids), nk)
        x <- extractClusters(sce, k, ks, verbose = FALSE)
        # check object naming
        expect_identical(length(x), nk)
        expect_identical(c(fsApply(x, identifier)), ks)
        # check no. of cells
        expect_identical(
            as.numeric(fsApply(x, nrow)), 
            as.numeric(table(kids)[ks]))
        # check expression matrices
        cs_by_k <- split(seq_len(ncol(sce)), kids)
        expect_true(all(vapply(ks, function(k) 
            all.equal(
                asinh(t(exprs(x[[k]]))/5), 
                assay(sce, "exprs")[, cs_by_k[[k]]]), 
            logical(1))))   
    })
})
