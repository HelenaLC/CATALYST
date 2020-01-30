context("gating")
library(flowCore)
library(SingleCellExperiment)

# construct SCE with some uniformly distributed data
n_cells <- 1e4
n_chs <- 20
chs <- paste0("ch", seq_len(n_chs))
x <- matrix(runif(n_cells * n_chs), n_cells, n_chs, dimnames = list(NULL, chs))
x <- fcs2sce(flowFrame(x), transform = FALSE)
assayNames(x) <- "exprs"

test_that("gateCytof() - type = 'rect'", {
    n <- 100 # number of cells that should be selected
    v <- 10  # expression value for cells to be selected
    i <- sample(n_chs, 2)   # sample 2 channels
    j <- sample(n_cells, n) # and 'n' cells
    # set these to v
    ij <- cbind(rep(i, n), rep(j, each = 2))
    y <- x; assay(y)[ij] <- v
    y <- gateCytof(y, chs[i], gate_id = "gate",
        type = "rect", xy = list(v-c(1,1), v+c(1,1)))
    expect_is(y, "SingleCellExperiment")
    expect_is(y$gate, "logical")
    expect_true(sum(y$gate) == n)
    expect_true(int_metadata(y)$gates$gate$type == "rect")
    # same ID with `overwrite = FALSE` should fail
    expect_error(gateCytof(y, chs[i], gate_id = "gate", 
        type = "quad", xy = c(1,1), overwrite = FALSE))
    # same ID with `overwrite = TRUE` should succeed
    expect_silent(gateCytof(y, chs[i], gate_id = "gate", 
        type = "quad", xy = c(1,1), overwrite = TRUE))
})

test_that("gateCytof() - type = 'elip'", {
    n <- 100 # number of cells that should be selected
    v <- 5  # expression value for cells to be selected
    i <- sample(n_chs, 2)   # sample 2 channels
    j <- sample(n_cells, n) # and 'n' cells
    # set these to v
    ij <- cbind(rep(i, n), rep(j, each = 2))
    y <- x; assay(y)[ij] <- v + runif(n)
    # should render message when 'k' is unspecified
    expect_message(gateCytof(y[, j], chs[i], type = "elip", xy = c(v,v)))
    y <- gateCytof(y, chs[i], type = "elip", k = 2, xy = c(v,v))
    expect_true(sum(y$gate1) == n)
})

test_that("gateCytof() - type = 'quad'", {
    n <- 100 # number of cells that should be selected
    v <- 10  # expression offset for cells to be selected
    # sample 2 channels and 'n' cells per quadrant
    i <- sample(n_chs, 2)    
    j <- sample(n_cells, 4*n)
    j <- split(j, seq_len(4))
    # shift cells by -/+10 for each quadrant
    s <- c(1, -1)
    s <- expand.grid(s, s)
    y <- x; for (k in seq_along(j)) {
        assay(y)[i[1], j[[k]]] <- s[k, 1] * v
        assay(y)[i[2], j[[k]]] <- s[k, 2] * v
    }
    y <- gateCytof(y, chs[i], gate_id = "gate",
        type = "quad", xy = c(-3, -3))
    expect_true(all(table(y$gate)[-1] == n))
})
