data(sample_ff, sample_key)
n_ids <- length(ids <- rownames(sample_key))
x <- prepData(sample_ff, by_time = FALSE)
x <- assignPrelim(x, sample_key, verbose = FALSE)

test_that("plotYields()", {
    expect_error(plotYields(x, which = "x"))
    # summary plot
    p <- plotYields(x, which = "0")
    expect_is(p, "ggplot")
    h <- p$layers[[1]]$data # histo data
    l <- p$layers[[2]]$data # lines data
    expect_true(nrow(h) == (n_cuts <- 101))
    expect_true(nrow(l) == n_cuts*n_ids)
    expect_true(all(table(l$bc_id) == n_cuts))
    expect_true(all(table(l$cutoff) == n_ids))
    expect_true(all(apply(h, 2, is.numeric)))
    expect_is(l$bc_id, "factor")
    expect_is(l$yield, "numeric")
    expect_is(l$cutoff, "numeric")
    expect_gte(min(l$yield), 0)
    expect_equal(min(l$yield), min(h$count))
    expect_equal(max(l$yield), max(h$count))
    # single ID
    p <- plotYields(x, which = "A1")
    expect_is(p, "ggplot")
    expect_true(all(apply(p$data, 2, is.numeric)))
    expect_gte(min(p$data$yield), 0)
    expect_lte(max(p$data$yield), 100)
    # multiple IDs
    n <- sample(seq_len(nrow(sample_key))[-1], 1)
    which <- sample(ids, n)
    m <- match(ids, which, nomatch = 0)
    p <- plotYields(x, which)
    expect_is(p, "list")
    expect_true(length(p) == n)
    expect_identical(names(p), which[m])
    # with separation cutoff estimates available
    y <- estCutoffs(x)
    which <- sample(ids, 3)
    p <- plotYields(y, which)
    cuts <- vapply(p[which], function(u) 
        as.numeric(u$layer[[4]]$data), numeric(1))
    expect_identical(cuts, metadata(y)$sep_cutoffs[which])
    # save output to .pdf
    foo <- plotYields(x, ids[1], 
        out_path = (dir <- tempdir()),
        out_name = (fn <- "foo"))
    expect_true(paste0(fn, ".pdf") %in% list.files(dir))
})

test_that("plotEvents()", {
    n <- 100 # number of events to plot
    expect_error(plotEvents(x, which = "x"))
    # single ID
    p <- plotEvents(x, which = "A1", n = n)
    expect_is(p, "ggplot")
    expect_true(nrow(p$data) == n * ncol(sample_key))
    expect_true(all(table(p$data$variable) == n))
    expect_equivalent(levels(p$data$variable), channels(x)[rowData(x)$is_bc])
    # multiple IDs
    n_in <- sample(seq_len(n_ids)[-1], 1)
    which <- sample(ids, n_in)
    m <- match(ids, which, nomatch = 0)
    p <- plotEvents(x, which, n = n)
    expect_is(p, "list")
    expect_true(all(sapply(p, is, "ggplot")))
    expect_true(length(p) == n_in)
    expect_identical(names(p), which[m])
    # all IDs
    p <- plotEvents(x, which = "all", n = n)
    expect_is(p, "list")
    expect_true(all(sapply(p, is, "ggplot")))
    expect_true(length(p) == n_ids)
    expect_identical(names(p), ids)
    # population(s) w/o events assigned
    ex_id <- sample(ids, (n_ex <- 1))
    in_ids <- sample(setdiff(ids, ex_id), (n_in <- 3))
    y <- x; y <- y[, y$bc_id != ex_id]
    expect_error(plotEvents(y, which = ex_id))
    p <- plotEvents(y, which = c(in_ids, ex_id), n = n)
    expect_is(p, "list")
    expect_true(all(sapply(p, is, "ggplot")))
    expect_true(length(p) == n_in)
    m <- match(ids, in_ids, nomatch = 0)
    expect_identical(names(p), in_ids[m])
    # save output to .pdf
    foo <- plotEvents(x, ids[1], n = n, 
        out_path = (dir <- tempdir()),
        out_name = (fn <- "foo"))
    expect_true(paste0(fn, ".pdf") %in% list.files(dir))
})

test_that("plotMahal()", {
    x <- applyCutoffs(estCutoffs(x))
    expect_error(plotMahal(x, which = "x"))
    expect_error(plotMahal(x, which = sample(ids, 2)))
    p <- plotMahal(x, which = "A1", n = 100)
    expect_is(p, "gtable")
})

test_that("plotSpillmat()", {
    data(ss_exp)
    bc_ms <- c(139, 141:156, 158:176)
    x <- prepData(ss_exp)
    x <- assignPrelim(x, bc_ms)
    x <- applyCutoffs(estCutoffs(x))
    x <- computeSpillmat(x)
    sm <- metadata(x)$spillover_matrix
    y <- x; metadata(y) <- list()
    p <- plotSpillmat(x)
    q <- plotSpillmat(y, sm)
    expect_is(p, "ggplot")
    expect_error(plotSpillmat(y))
    expect_identical(p$data, q$data)
})
