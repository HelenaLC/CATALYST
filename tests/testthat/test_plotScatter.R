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

test_that("plotScatter() - > 2 channels", {

})

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

#devtools::install('~/packages/CATALYST')
devtools::load_all('~/packages/CATALYST')

chs <- rownames(x)[1:2]
chs2 <- sample(rownames(x), 3)
type <- "rect"
group_by <- "bc_id"
xy <- list(c(4,4), c(8,8))
gate_id <- "x"
gate_id2 <- "y"
show_gate <- TRUE
show_perc <- TRUE
assay <- "exprs"
bins <- 100

y <- gateCytof(x0, chs, type = type, group_by = group_by, 
    overwrite = TRUE, xy = xy, gate_id = gate_id)
z <- gateCytof(y, rownames(x)[1:2], type = "rect", group_by = "bc_id", 
    overwrite = TRUE,xy = list(c(5,5), c(7,7)), gate_on = gate_id, gate_id = gate_id2)

gi <- int_metadata(z)$gates$y

plotScatter(z, chs = chs2, gate_id = gate_id,
    show_perc = TRUE, show_gate = TRUE, color_by = "selection")
plotScatter(z, chs = chs2, gate_id = gate_id2,
    show_perc = TRUE, show_gate = TRUE, color_by = "selection")

devtools::load_all()
plotScatter(z, chs2, gate_id = gate_id2, color_by = "selection")

plotScatter(y, chs = sample(rownames(x), 3), gate_id = "gate1", color_by = "selection")

plotScatter(x, chs = sample(rownames(x), 3))


devtools::load_all()
plotScatter(y, chs = sample(rownames(x), 3), show_perc = FALSE, gate_id = "gate1")

devtools::load_all()
gateCytof(x, rownames(x)[1:2], group_by = "bc_id", type = "rect", xy = list(c(1,1), c(4,4)))

devtools::load_all()
chs <- sample(rownames(x), 3)
plotScatter(z, chs = sample(rownames(x), 3), gate_id = "gate1")
plotScatter(z[, z$bc_id != "0"], chs = sample(rownames(x), 3), gate_id = "gate1")

plotScatter(z, chs = sample(rownames(x), 3), gate_id = "gate1", color_by = "selection")
plotScatter(z[, z$bc_id != "0"], chs = sample(rownames(x), 3), gate_id = "gate1", color_by = "selection")
