data(raw_data)
x <- prepData(raw_data)

test_that("normCytof()", {
    y <- normCytof(x, "dvs", k = 101, remove_beads = TRUE, overwrite = FALSE)
    sce <- "SingleCellExperiment"
    expect_is(y$data, sce)
    expect_is(y$beads, sce)
    expect_is(y$removed, sce)
    expect_is(y$scatter, "ggplot")
    expect_is(y$lines, "ggplot")
    expect_true(all(!y$data$remove))
    expect_true(all(!y$data$is_bead))
    expect_true(ncol(x) == ncol(y$data) + ncol(y$removed))

    z <- normCytof(x, "dvs", k = 101, remove_beads = FALSE, overwrite = FALSE)
    expect_identical(dim(z$data), dim(x))
    expect_identical(int_colData(z$data)$Time, int_colData(x)$Time)
    expect_equivalent(counts(z$data), counts(x))
    expect_lte(sum(z$data$is_bead), sum(z$data$remove))
    expect_true(sum(z$data$is_bead) == ncol(y$beads))
    expect_true(sum(z$data$is_bead | z$data$remove) == ncol(y$removed))
    
    # construct some mock data w/
    # signal descreasing over time
    set.seed(42)
    x <- prepData(raw_data)
    chs <- channels(x)
    bead_chs <- sample(chs, (n_beads <- 3))
    bead_ms <- .get_ms_from_chs(bead_chs)
    # amount time-drift (slope) & time points
    s <- -0.5 
    t <- seq(0, 10, l = (n <- 2e3))
    bl <- runif(n_beads, 5, 10) # baselines
    z <- outer(s*t, bl-s*max(t)/2, "+")
    # add time, DNA & some other channels
    z <- cbind(z, t, 0, 0, replicate((n_chs <- 5), runif(n)))
    # set random non-bead events
    i <- sample(n, (n_cells <- 800))
    j <- replicate(n_cells, sample(seq_len(n_beads), 1))
    z[cbind(i, j)] <- 0
    # set colnames to beads, dna, random channels
    dna <- c("Ir191Di", "Ir193Di")
    chs <- channels(x)
    chs <- sample(setdiff(chs, c(bead_chs, dna)), n_chs)
    colnames(z) <- c(bead_chs, "time", dna, chs) 
    # consruct SCE & apply normalization
    sce <- prepData(flowFrame(z))
    assay(sce, "exprs", FALSE) <- assay(sce)
    res <- normCytof(sce, bead_ms, k = 7, plot = FALSE, overwrite = TRUE)
    # check number of identified beads
    expect_equal(ncol(res$beads), n-n_cells)
    # fit LM model through normalized beads
    normed <- data.frame(t = t[setdiff(seq_len(n), i)],
        t(assay(res$beads, "counts")[bead_chs, ]))
    coefs <- vapply(bead_chs, function(u) {
        f <- as.formula(paste(u, "~ t"))
        coef(lm(f, data = normed))
    }, numeric(2))
    # LM fit intercepts should be similar to simulated baselines
    expect_true(all(abs(coefs[1, ]-bl) < 0.25))
    # LM fit slopes should be close to zero after normalization
    expect_true(all(abs(coefs[2, ]) < 0.1))
})

test_that("normCytof() - overwrite = TRUE", {
    i <- sample(ncol(x), 1e3)
    y <- normCytof(x[, i], beads = "dvs", k = 21, 
        overwrite = FALSE, transform = TRUE, plot = FALSE)
    z <- normCytof(x[, i], beads = "dvs", k = 21, 
        overwrite = TRUE, transform = TRUE, plot = FALSE)
    expect_true(!"ncounts" %in% assayNames(z$data))
    expect_true(2*length(assays(z$data)) == length(assays(y$data)))
    expect_identical(assay(y$data, "normcounts"), assay(z$data, "counts"))
    expect_identical(assay(y$data, "normexprs"), assay(z$data, "exprs"))
})

test_that("normCytof() - cofactor = NULL", {
    cfs <- sample(10, nrow(x), TRUE)
    names(cfs) <- channels(x)
    x <- prepData(raw_data[[1]], cofactor = cfs)
    i <- sample(ncol(x), 1e3)
    y <- normCytof(x[, i], beads = "dvs", k = 21, cofactor = NULL)
    expect_identical(int_metadata(y$data)$cofactor, cfs)
})
