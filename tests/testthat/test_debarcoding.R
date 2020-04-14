data(sample_ff, sample_key)
n_ids <- length(ids <- rownames(sample_key))
x <- prepData(sample_ff, by_time = FALSE)
x <- assignPrelim(x, sample_key, "counts", FALSE)

test_that("assignPrelim(); single-positive populations", {
    n_chs <- 10; n_cells <- 100
    x <- do.call(rbind, replicate(n_cells, diag(n_chs), simplify = FALSE))
    colnames(x) <- paste0("ch", seq_len(n_chs))
    x <- prepData(flowFrame(x), by_time = FALSE)
    x <- assignPrelim(x, seq_len(n_chs), "counts")
    expect_true(all(x$delta == 1))
    expect_true(all(table(x$bc_id) == n_cells))
    expect_equal(range(assay(x, "scaled")), c(0, 1))
})

test_that("assignPrelim(); 6-choose-3 debarcoding scheme", {
    expect_true(all(table(x$bc_id) == 1e3))
    expect_true(all(x$bc_id %in% rownames(sample_key)))
})

test_that("assignPrelim(); 
    debarcoding scheme with non-constant number of ones", {
    # simulate debarcoding scheme with 2-3/5 positive 
    bc_chs <- colnames(sample_key)
    n_bcs <- 5 # number of barcode channels
    n_pos <- seq_len(3)[-1] # number of possibly positive
    in_bcs <- sample(bc_chs, n_bcs) # sample barcode channels
    ex_bcs <- setdiff(bc_chs, in_bcs)
    key <- do.call(rbind, lapply(n_pos, function(n) {
        pos <- combn(n_bcs, n)
        n_ids <- choose(n_bcs, n)
        code <- t(replicate(n_ids, numeric(n_bcs)))
        i <- cbind(rep(seq_len(n_ids), each = n), c(pos))
        code[i] <- 1; return(code)
    }))
    ids <- paste0("id", seq_len(nrow(key)))
    rownames(key) <- ids
    colnames(key) <- in_bcs
    # sample new IDs & construct artificial data
    cs <- split(seq_len(n <- 1e3), sample(ids, n, replace = TRUE))
    ns <- vapply(cs, length, numeric(1))
    a <- lapply(ids, function(id) {
        pos <- key[id, ]
        neg <- rep(0, length(ex_bcs))
        replicate(ns[id], c(pos, neg))
    })
    m <- match(bc_chs, c(in_bcs, ex_bcs))
    a <- do.call(cbind, a)[m, ]
    y <- x[, seq_len(n)]
    assay(y, "exprs", withDimnames = FALSE) <- a
    # sample some IDs to be unassigned
    ex_ids <- sample(ids, 3)
    in_ids <- setdiff(ids, ex_ids)
    y <- assignPrelim(y, key[in_ids, ], verbose = FALSE)
    ns <- c("0" = sum(ns[ex_ids]), ns[in_ids])
    expect_equal(c(table(y$bc_id))[names(ns)], ns)
})

test_that("estCutoffs()", {
    x <- estCutoffs(x)
    y <- metadata(x)$sep_cutoffs
    expect_is(y, "numeric")
    expect_true(all(y[!is.na(y)] <= 1))
    expect_true(all(y[!is.na(y)] >= 0))
    expect_identical(names(y), rownames(sample_key))
    z <- x; z$bc_id <- NULL; expect_error(estCutoffs(z))
    z <- x; z$delta <- NULL; expect_error(estCutoffs(z))
})

test_that("applyCutoffs()", {
    y <- applyCutoffs(x, assay = "counts", mhl_cutoff = Inf, sep_cutoffs = 0)
    expect_identical(y$bc_id, x$bc_id)
    y <- applyCutoffs(x, assay = "counts", mhl_cutoff = 0, sep_cutoffs = 0)
    expect_true(all(y$bc_id == 0))
    y <- applyCutoffs(x, assay = "counts", sep_cutoffs = Inf)
    expect_true(all(y$bc_id == 0))
    y <- applyCutoffs(x, assay = "counts", mhl_cutoff = Inf, sep_cutoffs = 0.5)
    expect_equal(sum(y$bc_id == 0), sum(y$delta < 0.5))
    seps <- runif(nrow(sample_key), 0, 0.5); names(seps) <- ids
    y <- applyCutoffs(x, assay = "counts", mhl_cutoff = Inf, sep_cutoffs = seps)
    ds <- split(x$delta, x$bc_id)
    expect_equal(c(table(y$bc_id)[ids]), vapply(ids, 
        function(id) sum(ds[[id]] > seps[id]), numeric(1)))
})
