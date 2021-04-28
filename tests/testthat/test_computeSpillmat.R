library(flowCore)

# generate some dummy testdata with ss_exp as a template
data(ss_exp)
set.seed(12345)
n_cells <- 100

# generate a known spillover matrix
chs <- colnames(ss_exp)
sm <- diag(n_chs <- ncol(ss_exp))
colnames(sm) <- rownames(sm) <- chs
# add isotope, +/-1, +16 specific spillover
sm_ms <- .get_ms_from_chs(chs)
sm_mets <- .get_mets_from_chs(chs)

for (i in seq_along(sm_ms)) {
    ch_i <- chs[i]
    curmass <- sm_ms[i]
    curmet <- sm_mets[i]
    iso <- sm_mets == sm_ms[i]
    p1 <- sm_ms == sm_ms[i] + 1
    m1 <- sm_ms == sm_ms[i] - 1
    ox <- sm_ms == sm_ms[i] + 16
    spill_cols <- chs[(iso | p1 | m1 | ox) & (chs != chs[i])]
    sm[chs[i], spill_cols] <- runif(length(spill_cols)) / 20
}

# generate 100 cells that are positive for each channel
mat <- do.call(rbind, replicate(n_cells, diag(n_chs), simplify = FALSE))

# add noise & spillover
mat <- mat * (rnorm(n_chs * n_cells, 1, 0.01)) * 100
mat_sm <- mat %*% sm

# construct SCE & debarcode
sce <- prepData(flowFrame(mat_sm), by_time = FALSE)
sce <- assignPrelim(sce, sm_ms, verbose = FALSE)

# set barcode IDs to constructed IDs
sce$bc_id <- rep(as.character(sm_ms), n_cells)
    
for (i in c("default", "all")) 
    for (m in c("default", "classic")) 
test_that(sprintf("interactions = '%s', method = '%s'", i, m), {
    # estimate spillover matrix
    sce <- computeSpillmat(sce, interactions = i, method = m)
    sm_est <- metadata(sce)$spillover_matrix
    expect_equal(sm_est, sm, info = c(
        "test if spillover matrix can be reconstructed",
        "from simulated single cells without noise."))
    
    # repeat with a constant amount of background
    sce <- prepData(flowFrame(mat_sm + 20), by_time = FALSE)
    sce <- assignPrelim(sce, sm_ms, verbose = FALSE)
    sce$bc_id <- rep(as.character(sm_ms), n_cells)
    sce <- computeSpillmat(sce)
    sm_est <- metadata(sce)$spillover_matrix
    expect_equal(sm_est, sm, info = c(
        "test if spillover matrix can be reconstructed",
        "from simulated single cells without noise but with background."))
    
    noise <- rnorm(n_chs * n_cells, 5, 1)
    sce <- prepData(flowFrame(mat_sm * noise), by_time = FALSE)
    sce <- assignPrelim(sce, sm_ms, verbose = FALSE)
    sce$bc_id <- rep(as.character(sm_ms), n_cells)
    sce <- computeSpillmat(sce)
    sm_est <- metadata(sce)$spillover_matrix
    expect_equal(sm_est, sm, info = c(
        "test if spillover matrix can be reconstructed",
        "from simulated single cells without noise but with background."))
})
