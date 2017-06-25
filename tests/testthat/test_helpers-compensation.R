# ==============================================================================
# unit tests for compensation helpers
# ------------------------------------------------------------------------------
test_that("make_symetric() works flawlessly.", {
    dims <- sample(1:10, 2)
    min <- which.min(dims)
    max <- which.max(dims)
    nms <- as.character(runif(dims[max]))
    dimNms <- vector("list", 2)
    dimNms[[min]] <- sample(nms, dims[min])
    dimNms[[max]] <- nms
    m <- matrix(seq_len(dims[1]*dims[2]), dims[1], dims[2], dimnames=dimNms)
    m <- make_symetric(m)
    
    expect_true(!any(is.na(m)))
    expect_true(is.matrix(m))
    expect_true(is.numeric(m))
    expect_true(length(unique(dim(m))) == 1)
    expect_true(all.equal(rownames(m), colnames(m)))
    expect_error(make_symetric("NaN"))
})

test_that("get_spill_cols() works impeccably", {
    n <- 25
    valid_cols <- seq_len(n)
    ms <- sample(139:176, n)
    mets <- sample(size=n, x=c(
        "Ba","La","Ce","Pr","Nd","Nd","Nd","Nd","Nd","Sm","Nd","Sm","Nd",
        "Eu","Sm","Eu","Sm","Gd","Gd","Gd","Gd","Tb","Gd","Dy","Dy","Dy",
        "Dy","Ho","Er","Er","Er","Tm","Er","Yb","Yb","Yb","Yb","Lu","Yb"))
    l <- get_spill_cols(ms, mets)
    
    expect_true(is.list(l))
    expect_true(all(unlist(lapply(l, function(i) i %in% valid_cols))))
    expect_true(all(sapply(l, function(i) length(unique(i)) == length(i))))
})
    