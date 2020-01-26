context("compensation")

data(ss_exp)
library(flowCore)
library(SingleCellExperiment)

test_that("compCytof()", {
    # test if artificially added spillover can be removed
    # generate a dummy spillover matrix
    n_chs <- ncol(ss_exp)
    sm <- diag(1, n_chs, n_chs)
    inds <- cbind(seq_len(n_chs - 1), seq_len(n_chs)[-1])
    sm[inds] <- seq(0.01, 0.3, l = n_chs - 1)
    rownames(sm) <- colnames(sm) <- colnames(ss_exp)
    
    # add artificial spillover
    ss_spill <- ss_exp
    exprs(ss_spill) <- exprs(ss_spill) %*% sm
    # construct SCE & compensate
    sce <- fcs2sce(ss_spill, by_time = FALSE)
    sce <- compCytof(sce, sm, method = "flow")
    s <- adaptSpillmat(sm, rownames(sce)) 
    # get raw & compensated data
    ref <- t(exprs(ss_exp))
    comped <- assay(sce, "counts_comped")
    expect_equal(ref, comped, tolerance = 10^-10, 
        info = "compCytof() corrects artificial spillover.")
    
    # to test if the spillover addition was correct, 
    # the data gets also compensated using 'flowCore' directly
    comped_flowCore <- t(exprs(compensate(ss_spill, sm)))
    expect_equal(ref, comped_flowCore, tolerance = 10^-10, 
        info = "spillover is simulated correctly.")
    
    # test nnls
    sce <- compCytof(sce, sm, method = "nnls")
    comped_nnls <- assay(sce, "counts_comped")
    expect_equal(ref, comped_nnls, tolerance = 10^-10, 
        info = "compCytof() NNLS-method corrects artificial spillover.")
})
