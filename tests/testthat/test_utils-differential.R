suppressPackageStartupMessages({
    library(diffcyt)
    library(flowCore)
    library(SingleCellExperiment)
})

x <- .toySCE()
es <- assay(x, "exprs")

test_that(".get_features()", {
    expect_error(.get_features(x, "x"))
    expect_error(.get_features(x, c("x", rownames(x))))
    expect_identical(.get_features(x, NULL), rownames(x))
    valid <- c("type", "state", "none")
    classes <- sample(valid, nrow(x), TRUE)
    rowData(x)$marker_class <- classes
    for (c in valid) 
        expect_equivalent(
            .get_features(x, c), 
            rownames(x)[classes == c])
    c <- sample(valid, 1)
    rowData(x)$marker_class <- c
    expect_error(.get_features(x, setdiff(valid, c)))
})

test_that(".get_shapes()", {
    for (cd_var in names(colData(x))) {
        shapes <- .get_shapes(x, cd_var)
        expect_true(length(shapes) == nlevels(x[[cd_var]]))
    }
    y <- x; y$foo <- factor(seq_len(100))
    expect_message(shapes <- .get_shapes(y, "foo"))
    expect_true(is.null(shapes))
})

test_that(".check_pal()", {  
    expect_error(.check_pal("blue")) 
    expect_error(.check_pal(letters))
    expect_silent(.check_pal(colors()))
    expect_error(.check_pal(c("x", colors())))
})

test_that("guessPanel()", {
    expect_error(guessPanel("x"))
    ff <- get(data("sample_ff"))
    ps0 <- parameters(ff)
    # channel & antigen in 'name' & 'desc'
    expect_is(df <- guessPanel(ff), "data.frame")
    expect_true(nrow(df) == ncol(ff))
    expect_true(all(df$fcs_colname == ps0@data$name))
    expect_true(all(df$antigen == ps0@data$desc))
    expect_true(all(df$desc == ps0@data$desc))
    # 'desc' of the from channel_antigen
    ps <- ps0
    ps@data$desc <- with(ps@data, paste(name, desc, sep = "_"))
    parameters(ff) <- ps
    expect_is(df <- guessPanel(ff), "data.frame")
    expect_true(all(df$fcs_colname == ps@data$name))
    expect_true(all(df$antigen == ps0@data$desc))
    expect_true(all(df$desc == ps0@data$name))
})

test_that(".split_cells() by 1 factor", {
    for (by in c("sample_id", "cluster_id")) {
        cs <- .split_cells(x, by)
        expect_identical(names(cs), levels(x[[by]]))
        expect_identical(
            as.numeric(vapply(cs, length, numeric(1))),
            as.numeric(table(x[[by]])))
    }
})

test_that(".split_cells() by 2 factors", {
    for (by in list(
        c("sample_id", "cluster_id"), 
        c("cluster_id", "sample_id"))) {
        cs <- .split_cells(x, by)
        expect_identical(names(cs), levels(x[[by[1]]]))
        expect_true(all(vapply(cs, function(u) 
            identical(names(u), levels(x[[by[2]]])), 
            logical(1))))
    }
})

test_that(".agg() by 1 factor", {
    for (by in c("sample_id", "cluster_id")) {
        for (fun in c("sum", "mean", "median")) {
            pb <- .agg(x, by, fun)
            replicate(10, {
                i <- sample(seq_len(nrow(x)), 1)
                j <- sample(levels(x[[by]]), 1)
                expect_equal(pb[i, j], get(fun)(es[i, x[[by]] == j]))
            })
        }
    }
})

test_that(".agg() by 2 factors", {
    for (by in list(
        c("sample_id", "cluster_id"), 
        c("cluster_id", "sample_id"))) {
        for (fun in c("sum", "mean", "median")) {
            pb <- .agg(x, by, fun, "exprs")
            expect_identical(names(pb), levels(x[[by[1]]]))
            expect_true(all(vapply(pb, function(u) 
                identical(colnames(u), levels(x[[by[2]]])), 
                logical(1))))
            replicate(10, {
                a <- sample(levels(x[[by[1]]]), 1)
                b <- sample(levels(x[[by[2]]]), 1)
                i <- sample(nrow(x), 1)
                j <- x[[by[1]]] == a & x[[by[2]]] == b
                expect_equal(pb[[a]][i, b], get(fun)(es[i, j]))
            })
        }
    }
})

.diffcyt <- function(type = c("DA", "DS"), n = 20, m = 10, ns = 4) {
    data(PBMC_fs, PBMC_panel, PBMC_md)
    x <- prepData(PBMC_fs, PBMC_panel, PBMC_md)
    
    ps <- sample(levels(x$patient_id), 2)
    x <- filterSCE(x, patient_id != ps)
    
    design <- createDesignMatrix(ei(x), cols_design = 2)
    contrast <- createContrast(c(0, 1))
    
    kids <- sample(c("a", "b"), ncol(x), TRUE)
    x$cluster_id <- factor(kids)
    metadata(x)$cluster_codes <- "foo"
    
    l <- switch(type, 
        DA = diffcyt(x, 
            design = design, contrast = contrast, verbose = FALSE,
            analysis_type = "DA", method_DA = "diffcyt-DA-edgeR"),
        DS = diffcyt(x, 
            design = design, contrast = contrast, verbose = FALSE,
            analysis_type = "DS", method_DS = "diffcyt-DS-limma"))
    
    return(rowData(l$res))
}
