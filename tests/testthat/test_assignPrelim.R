test_that("assignPrelim() works faultlessly.", {
    data <- matrix(0, nrow=100, ncol=10, 
        dimnames=list(NULL, paste0("channel", 1:10)))
    inds <- sample(1:10, 100, replace=TRUE)
    data[cbind(1:100, inds)] <- 1
    ff <- flowFrame(data)
    ids <- sample(1:10, 8)
    re <- assignPrelim(ff, ids)
    
    expect_is(re, "dbFrame")
    expect_true(all(rownames(bc_key(re)) == ids))
    expect_true(all(unique(bc_ids(re)) %in% c(0, ids)))
    expect_true(sum(inds == ids[1]) == sum(bc_ids(re) == ids[1]))
    expect_true(all.equal(dim(normed_bcs(re)), c(100, 8)))
    expect_true(sum(inds %in% ids) == sum(bc_ids(re) != 0))
    expect_error(assignPrelim(ff, 0))
    expect_error(assignPrelim("Not FCS", 1:10))
    expect_true(all.equal(c(exprs(ff)), c(exprs(re))))
})