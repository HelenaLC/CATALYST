test_that("compCytof() works.", {
    # Tests if artificail added spillover can be removed
    data(ss_exp)
    
    # generate a dummy spillover matrix
    ncol = flowCore::ncol(ss_exp)
    sm = diag(1, ncol, ncol)
    sm[cbind(1:(ncol-1),2:(ncol))] = 0.1
    colnames(sm) <- rownames(sm) <- flowCore::colnames(ss_exp)
    
    # add the spillover
    ss_morespill = ss_exp
    flowCore::exprs(ss_morespill)  = flowCore::exprs(ss_exp) %*% (sm)
    # compensate
    ss_morespill_comp = CATALYST::compCytof(ss_morespill, sm, method='flow')
    
    s = adaptCompensationSpillmat(sm, flowCore::colnames(ss_morespill)) 
    # check the data
    dat_orig = flowCore::exprs(ss_exp)
    dat_spill =  flowCore::exprs(ss_morespill)
    dat_corr = flowCore::exprs(ss_morespill_comp)
    # to test if the spillover addition was correct, the data gets also compensated using
    # flowcore directly
    dat_corr_flowcore = flowCore::exprs(flowCore::compensate(ss_morespill,sm))
    
    testthat::expect_equal(dat_orig, dat_corr_flowcore, tolerance=10^-10, info = 'Tests if the spillover was correctly simulated.')
    
    testthat::expect_equal(dat_orig, dat_corr, tolerance=10^-10, info = 'Tests if the compCytof corrects artificial spillover.')
    
    # test nnls
    dat_corr_nnls = flowCore::exprs(CATALYST::compCytof(ss_morespill, sm, method='nnls'))
    testthat::expect_equal(dat_orig, dat_corr_nnls, tolerance=10^-10, info = 'Tests if the compCytof nnls corrects artificial spillover.')
})