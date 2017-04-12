# read input FCS files, render box #2
observe({
    x <- input$fcsNorm
    if (is.null(x)) return()
    for (i in seq_len(nrow(x)))
        vals$ffsNorm[[i]] <- flowCore::read.FCS(x[[i, "datapath"]])
    output$box2 <- renderUI(box2)
    js$collapse("box_1")
})

# checkbox "Normalize to median level of current files"
observe({
    x <- input$box_NormToCurrent
    if (is.null(x) || x == 0) 
        return()
    updateCheckboxInput(session, "box_uploadNormTo", value=FALSE)
    if (length(vals$ffsNorm) == 1) {
        vals$ffsNormTo <- vals$ffsNorm[[1]]
    } else {
        vals$ffsNormTo <- CATALYST::concatFCS(vals$ffsNorm)
    }
    js$collapse("box_2")
})

# checkbox "Upload FCS file(s) of beads to normalize to"
observe({
    x <- input$box_uploadNormTo
    if (is.null(x) || x == 0) 
        return()
    updateCheckboxInput(session, "box_NormToCurrent", value=FALSE)
})

output$input_NormTo <- renderUI({
    x <- input$box_uploadNormTo
    if (is.null(x) || x == 0)
        return()
    fileInput("fcsNormTo", NULL, accept=".fcs")
})

# read FCS of beads to normalize to
observe({
    x <- input$fcsNormTo
    if (is.null(x)) return()
    if (nrow(x) > 1) {
        vals$ffsNormTo <- vector(mode="list", length=nrow(x))
        for (i in seq_len(nrow(x)))
            vals$ffsNormTo[[i]] <- flowCore::read.FCS(x[[i, "datapath"]])
        vals$ffsNormTo <- concatFCS(vals$ffsNormTo)
    } else {
        vals$ffsNormTo <- flowCore::read.FCS(input$fcsNormTo$datapath)
    }

    js$collapse("box_2")
})

# initialize logical vector of bead indices for each flowFrame
# and list for normalized flowFrame(s), render box #3
observe({
    if (is.null(vals$ffsNormTo)) 
        return()
    n1 <- length(vals$ffsNorm)
    if (input$box_NormToCurrent == 1) {
        vals$beadInds  <- vector(mode="list", length=n1)
    } else {
        n2 <- length(vals$ffsNormTo)
        vals$beadInds  <- vector(mode="list", length=sum(n1, n2))
    }
    vals$ffsNormed <- vector(mode="list", length=n1)
    output$box3 <- renderUI(box3)
})

# keep track of currently selected sample
observe({
    x <- input$select_sample
    if (is.null(x)) return()
    vals$smpl <- which(c(input$fcsNorm$name, input$fcsNormTo$name) == x)
})

# bead gating
observe({
    if (is.null(input$gate1)
        || is.null(input$gate2)
        || is.null(input$gate3)
        || is.null(input$gate4)
        || is.null(input$gate5))
        return()
    
    es <- flowCore::exprs(c(vals$ffsNorm, vals$ffsNormTo)[[vals$smpl]])
    chs <- colnames(es)
    
    # ----- get brushed points -----
    gate1 <- brushedPoints(brush=input$gate1, allRows=TRUE,
        df=data.frame(asinh(es/input$input_cfGating)),
        xvar=chs[vals$beadCols[1]], yvar=chs[vals$dnaCols[1]])
    gate2 <- brushedPoints(brush=input$gate2, allRows=TRUE,
        df=data.frame(asinh(es/input$input_cfGating)),
        xvar=chs[vals$beadCols[2]], yvar=chs[vals$dnaCols[1]])
    gate3 <- brushedPoints(brush=input$gate3, allRows=TRUE,
        df=data.frame(asinh(es/input$input_cfGating)),
        xvar=chs[vals$beadCols[3]], yvar=chs[vals$dnaCols[1]])
    gate4 <- brushedPoints(brush=input$gate4, allRows=TRUE,
        df=data.frame(asinh(es/input$input_cfGating)),
        xvar=chs[vals$beadCols[4]], yvar=chs[vals$dnaCols[1]])
    gate5 <- brushedPoints(brush=input$gate5, allRows=TRUE,
        df=data.frame(asinh(es/input$input_cfGating)),
        xvar=chs[vals$beadCols[5]], yvar=chs[vals$dnaCols[1]])
    
    gate1 <- gate1[, ncol(gate1)]
    gate2 <- gate2[, ncol(gate2)]
    gate3 <- gate3[, ncol(gate3)]
    gate4 <- gate4[, ncol(gate4)]
    gate5 <- gate5[, ncol(gate5)]
    # -----
    
    # geat indices of events falling in intersection of all gates
    vals$beadInds[[vals$smpl]] <- gate1 & gate2 & gate3 & gate4 & gate5
})

# gating yield valueBox
observe({
    y <- vals$beadInds[[vals$smpl]]
    if (is.null(y)) { 
        yield <- "0.00"
        smplNm <- "of events gated"
    } else {
        x <- c(vals$ffsNorm, vals$ffsNormTo)[[vals$smpl]]
        yield <- sprintf("%2.2f", sum(y)/nrow(x)*100)
        smplNm <- c(input$fcsNorm$name, input$fcsNormTo$name)[vals$smpl]
    }
    output$text_gatingYield <- renderText(yield)
    output$text_smplNm <- renderText(smplNm)
})

# normalization
observe({
    # check that all samples have been gated
    if (any(vapply(vals$beadInds, is.null, logical(1))))
        return()
    
    n1 <- length(vals$ffsNorm)
    n2 <- length(vals$ffsNormTo)
    if (input$box_NormToCurrent == 1) {
        inds <- vals$beadInds[1:n1]
    } else {
        inds <- vals$beadInds[(n1+1):n2]
    }
    
    # get baseline values from ffsNormTo and normalize ffsNorm
    bl <- getBaseline(vals$ffsNormTo, vals$beads, unlist(inds))
    for (i in seq_along(vals$ffsNorm)) {
        vals$ffsNormed[[i]] <- normCytof(
            vals$ffsNorm[[i]], vals$beads, vals$beadInds[[i]], bl
        )
    }
    output$box_smoothedBeads <- renderUI(box_smoothedBeads)
})

observe({
    if (any(vapply(vals$ffsNormed, is.null, logical(1))))
        return()
    if (vals$smpl > length(vals$ffsNorm)) 
        return()
    if (sum(vals$beadInds[[vals$smpl]]) == 0) {
        showNotification("No beads have been gated for this sample 
            Please consider re-gating.", type="error", closeButton=FALSE)
        return()
    }
    es <- flowCore::exprs(vals$ffsNorm[[vals$smpl]])
    es <- es[vals$beadInds[[vals$smpl]], ]
    
    normedEs <- flowCore::exprs(vals$ffsNormed[[vals$smpl]])
    normedEs <- normedEs[vals$beadInds[[vals$smpl]], ]
    
    smoothedBeads <- data.frame(
        es[, vals$timeCol],
        vapply(vals$beadCols, function(i)
            stats::runmed(es[, i], 501, "constant"),
        numeric(sum(vals$beadInds[[vals$smpl]]))))
    
    smoothedBeadsNormed <- data.frame(
        es[, vals$timeCol],
        vapply(vals$beadCols, function(i)
            stats::runmed(normedEs[, i], 501, "constant"),
        numeric(sum(vals$beadInds[[vals$smpl]]))))
    
    output$plot_smoothedBeads <- renderPlot(
        outPlots2(smoothedBeads, smoothedBeadsNormed, NULL))
})

# ··············································································
# bead vs. dna scatters
# ··············································································
observe({
    if (is.null(vals$beadCols) || is.null(vals$dnaCols)) 
        return()
    (input$button_viewGating)
    es <- flowCore::exprs(c(vals$ffsNorm, vals$ffsNormTo)[[vals$smpl]])

    output$plot_beadScatter1 <- renderPlot({ 
        plotScatter(
            es=es,
            x=vals$beadCols[1], 
            y=vals$dnaCols[1],
            cf=isolate(input$input_cfGating), 
            n=isolate(input$input_nGating))
    })
    output$plot_beadScatter2 <- renderPlot({ 
        plotScatter(
            es=es,
            x=vals$beadCols[2], 
            y=vals$dnaCols[1],
            cf=isolate(input$input_cfGating), 
            n=isolate(input$input_nGating))
    })
    output$plot_beadScatter3 <- renderPlot({ 
        plotScatter(
            es=es,
            x=vals$beadCols[3], 
            y=vals$dnaCols[1],
            cf=isolate(input$input_cfGating), 
            n=isolate(input$input_nGating))
    })
    output$plot_beadScatter4 <- renderPlot({ 
        plotScatter(
            es=es,
            x=vals$beadCols[4], 
            y=vals$dnaCols[1],
            cf=isolate(input$input_cfGating), 
            n=isolate(input$input_nGating))
    })
    output$plot_beadScatter5 <- renderPlot({ 
        plotScatter(
            es=es,
            x=vals$beadCols[5], 
            y=vals$dnaCols[1],
            cf=isolate(input$input_cfGating), 
            n=isolate(input$input_nGating))
    })
})