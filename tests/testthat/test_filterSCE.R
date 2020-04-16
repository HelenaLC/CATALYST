sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
sce <- cluster(sce, verbose = FALSE)

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))

test_that("filterSCE", {
    k <- sample(kids, sample(nk, 1))
    y <- filterSCE(sce, cluster_id %in% k)
    expect_identical(ncol(y), sum(sce$cluster_id %in% k))
    
    s <- sample(sids, sample(ns, 1))
    y <- filterSCE(sce, sample_id %in% s)
    expect_identical(ncol(y), sum(sce$sample_id %in% s))
    expect_identical(ei(y), droplevels(filter(ei(sce), sample_id %in% s)))
    expect_equal(sum(ei(y)$n_cells), ncol(y))
    m <- match(levels(y$sample_id), ei(y)$sample_id)
    expect_equal(ei(y)$n_cells[m], as.numeric(table(y$sample_id)))
    
    kid <- sample(names(cluster_codes(sce))[-1], 1)
    kids <- cluster_ids(sce, k = kid)
    k <- sample(levels(kids), sample(nlevels(kids), 1))
    s <- sample(sids, sample(ns, 1))
    y <- filterSCE(sce, k = kid, cluster_id %in% k, sample_id %in% s)
    expect_identical(ncol(y), sum(sce$sample_id %in% s & kids %in% k))
    expect_equal(sum(ei(y)$n_cells), ncol(y))
    m <- match(ei(y)$sample_id, levels(y$sample_id))
    expect_equal(ei(y)$n_cells, as.numeric(table(y$sample_id)[m]))
})
