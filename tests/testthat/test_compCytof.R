test_that("compCytof() is flowSOM", {
    # tests if artificially added spillover can be removed
    data(ss_exp)
    
    # generate a dummy spillover matrix
    ncol <- flowCore::ncol(ss_exp)
    sm <- diag(1, ncol, ncol)
    inds <- cbind(1:(ncol-1), 2:(ncol))
    sm[inds] <- seq(0.01, 0.3, length.out=ncol-1)
    colnames(sm) <- rownames(sm) <- flowCore::colnames(ss_exp)
    
    # add spillover
    ss_morespill <- ss_exp
    flowCore::exprs(ss_morespill) <- flowCore::exprs(ss_exp) %*% sm
    # compensate
    ss_morespill_comp <- CATALYST::compCytof(ss_morespill, sm, method="flow")
    
    s <- adaptSpillmat(sm, flowCore::colnames(ss_morespill)) 
    # check the data
    dat_orig <- flowCore::exprs(ss_exp)
    dat_spill <- flowCore::exprs(ss_morespill)
    dat_corr <- flowCore::exprs(ss_morespill_comp)
    
    # to test if the spillover addition was correct, 
    # the data gets also compensated using flowCore directly
    dat_corr_flowCore <- flowCore::compensate(ss_morespill, sm)
    dat_corr_flowCore <- flowCore::exprs(dat_corr_flowCore)
    
    testthat::expect_equal(dat_orig, dat_corr_flowCore, tolerance=10^-10, 
        info="Tests if the spillover was correctly simulated.")
    testthat::expect_equal(dat_orig, dat_corr, tolerance=10^-10, 
        info="Tests if the compCytof corrects artificial spillover.")
    
    # test nnls
    dat_corr_nnls <- CATALYST::compCytof(ss_morespill, sm, method='nnls')
    dat_corr_nnls <- flowCore::exprs(dat_corr_nnls)
    testthat::expect_equal(dat_orig, dat_corr_nnls, tolerance=10^-10, 
        info="Tests if the compCytof nnls corrects artificial spillover.")
})
