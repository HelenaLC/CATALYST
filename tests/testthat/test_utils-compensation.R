test_that(".check_sm()", {
    ns <- vapply((l <- isotope_list), length, numeric(1))
    n <- length(chs <- paste0(rep.int(names(l), ns), unlist(l)))
    x <- matrix(runif(n*n), n, n, dimnames = list(chs, chs))
    diag(x) <- 1; expect_silent(.check_sm(x))
    # diagonal entry != 1
    y <- x; y[1, 1] <- 0; expect_error(.check_sm(y))
    # entries out of range
    y <- x; y[1, 2] <- -1; expect_error(.check_sm(y))
    y <- x; y[1, 2] <- 10; expect_error(.check_sm(y))
    # spilling != receiving channels
    y <- x; rownames(y)[1] <- "Abc100"; expect_error(.check_sm(y))
    # invalid metal-mass combination
    colnames(y)[1] <- "Abc100"; expect_error(.check_sm(y))
    # modifying isotope list accordingly
    # should make SM valid again
    l <- c(l, Abc = 100); expect_silent(.check_sm(y, l))
})

test_that(".make_symetric()", {
    dims <- sample(seq_len(10), 2)
    min <- which.min(dims)
    max <- which.max(dims)
    nms <- as.character(runif(dims[max]))
    dimNms <- vector("list", 2)
    dimNms[[min]] <- sample(nms, dims[min])
    dimNms[[max]] <- nms
    m <- matrix(seq_len(dims[1]*dims[2]), dims[1], dims[2], dimnames=dimNms)
    m <- .make_symetric(m)
    
    expect_true(!any(is.na(m)))
    expect_true(is.matrix(m))
    expect_true(is.numeric(m))
    expect_true(length(unique(dim(m))) == 1)
    expect_true(all.equal(rownames(m), colnames(m)))
})

test_that("..get_spill_chs()", {
    n <- 25
    valid_cols <- seq_len(n)
    ms <- sample(seq(139, 176), n)
    mets <- sample(size=n, x=c(
        "Ba","La","Ce","Pr","Nd","Nd","Nd","Nd","Nd","Sm","Nd","Sm","Nd",
        "Eu","Sm","Eu","Sm","Gd","Gd","Gd","Gd","Tb","Gd","Dy","Dy","Dy",
        "Dy","Ho","Er","Er","Er","Tm","Er","Yb","Yb","Yb","Yb","Lu","Yb"))
    l <- .get_spill_chs(ms, mets)
    
    expect_true(is.list(l))
    expect_true(all(unlist(lapply(l, function(i) i %in% valid_cols))))
    expect_true(all(vapply(l, function(i) length(unique(i)) == length(i), logical(1))))
})

compare_sm_nms <- function(sm, nms){
    expect_true(all(colnames(sm) == nms))
    expect_true(all(rownames(sm) == nms))
}

test_that("adaptSpillmat()", {
    data(ss_exp)
    # generate a dummy spillover matrix
    ncol <- flowCore::ncol(ss_exp)
    sm <- diag(1, ncol, ncol)
    inds <- cbind(seq_len(ncol-1), seq_len(ncol)[-1])
    sm[inds] <- seq(0.01, 0.3, length.out = ncol-1)
    colnames(sm) <- rownames(sm) <- flowCore::colnames(ss_exp)
    chs <- flowCore::colnames(ss_exp)
    
    # case a subset of the channels
    inds <- seq_len(10)
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
