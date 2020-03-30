context("gating")
data(sample_ff, sample_key)
x <- fcs2sce(sample_ff)
# sample only a couple IDs for testing
ids <- sample(rownames(sample_key), 3)
x <- assignPrelim(x, sample_key[ids, ], verbose = FALSE)

# subset 2 channels & 'ids'
chs <- sample(rownames(x), 2)
x <- x[, x$bc_id %in% ids]
es <- t(assay(x[chs, ], "exprs"))
es <- data.frame(es, check.names = FALSE)

args <- list(
    list(type = "rect", geom = "GeomRect", xy = list(c(5,6), c(7,7))),
    list(type = "elip", geom = "GeomPath", xy = c(6,6)))
test_that("plotScatter() - gates with grouping", {
    for (i in args) {
        # apply gate
        y <- gateCytof(x, chs, group_by = "bc_id",
            k = 1, type = i$type, xy = i$xy)
        # scatter without gate
        p <- plotScatter(y, chs)
        expect_is(p, "ggplot")
        expect_identical(p$data[, chs], es)
        expect_true(length(p$facet$params) == 0)
        # scatter with facetting but without gate
        p <- plotScatter(y, chs, gate_id = "gate1", show_gate = FALSE)
        expect_is(p, "ggplot")
        expect_true(names(p$facet$params$facets) == "bc_id")
        ns <- vapply(split(p$data, p$data$bc_id), nrow, numeric(1))
        expect_equivalent(ns, c(table(x$bc_id)))
        # scatter with facetting & gate
        p <- plotScatter(y, chs, gate_id = "gate1", show_gate = TRUE)
        expect_is(p, "ggplot")
        expect_true(names(p$facet$params$facets) == "bc_id")
        ns <- vapply(split(p$data, p$data$bc_id), nrow, numeric(1))
        expect_equivalent(ns, c(table(x$bc_id)))
        expect_is(p$layers[[2]]$geom, i$geom)
        expect_identical(p$layers[[2]]$data, int_metadata(y)$gates$gate$data)
    }
})
