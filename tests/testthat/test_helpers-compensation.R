# ==============================================================================
# unit tests for compensation helpers
# ------------------------------------------------------------------------------
test_that("make_symetric() works flawlessly", {
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


compare_sm_nms <- function(sm, nms){
    expect_true(all(colnames(sm) == nms))
    expect_true(all(rownames(sm) == nms))
}

test_that("adaptSpillmat() is magnificent", {
    data(ss_exp)
    # generate a dummy spillover matrix
    ncol <- flowCore::ncol(ss_exp)
    sm <- diag(1, ncol, ncol)
    inds <- cbind(1:(ncol-1),2:(ncol))
    sm[inds] <- seq(0.01, 0.3, length.out=ncol-1)
    colnames(sm) <- rownames(sm) <- flowCore::colnames(ss_exp)
    chs <- flowCore::colnames(ss_exp)
    
    # case a subset of the channels
    inds <- 1:10
    sm_ad <- adaptSpillmat(sm, chs[inds])
    compare_sm_nms(sm_ad, chs[inds])
    # test that underlying spillover structure is not altered
    expect_equal(sm[chs[inds], chs[inds]], sm_ad[chs[inds], chs[inds]])
    
    # case a new channel
    chs2 <- c(chs, "Ce141Di")
    sm_ad <- adaptSpillmat(sm, chs2)
    compare_sm_nms(sm_ad, chs2)
    # test that underlying spillover structure is not altered
    expect_equal(sm[chs, chs], sm_ad[chs, chs])
    
    # case a multiple new channels
    chs2 <- c(chs, "Filename", "Ce141Di", "Time", "Ce139i")
    sm_ad <- adaptSpillmat(sm, chs2)
    compare_sm_nms(sm_ad, chs2)
    # test that underlying spillover structure is not altered
    expect_equal(sm[chs, chs], sm_ad[chs, chs])
    
    # case duplicated masses
    chs2 <- c(chs, "Ce141Di", "Pr140Di")
    sm_ad <- adaptSpillmat(sm, chs2)
    compare_sm_nms(sm_ad, chs2)
    # Pr140Di should now receive spillover as Ce140Di
    expect_equal(sum(sm_ad[, "Ce140Di"]), sum(sm_ad[, "Pr140Di"]))
    # Pr140 should not emmit any spillover
    expect_equal(sum(sm_ad["Pr140Di", ]),1)
    # test that underlying spillover structure is not altered
    expect_equal(sm[chs, chs], sm_ad[chs, chs])
    
    # case random channel order
    chs3 <- sample(chs2)
    sm_ad <- adaptSpillmat(sm, chs3)
    compare_sm_nms(sm_ad, chs3)
    # Pr140Di should now receive spillover as Ce140Di
    expect_equal(sum(sm_ad[, "Ce140Di"]), sum(sm_ad[, "Pr140Di"]))
    # Pr140 should not emmit any spillover
    expect_equal(sum(sm_ad["Pr140Di", ]),1)
    # test that underlying spillover structure is not altered
    expect_equal(sm[chs, chs], sm_ad[chs, chs])
    
    # case duplicated masses, random channel order and subset
    chs4 <- sample(chs3, 10)
    sm_ad <- adaptSpillmat(sm, chs4)
    compare_sm_nms(sm_ad, chs4)
    # test that underlying spillover structure is not altered
    chs4_orig <- chs[chs %in% chs4]
    expect_equal(sm[chs4_orig, chs4_orig], sm_ad[chs4_orig, chs4_orig])
})