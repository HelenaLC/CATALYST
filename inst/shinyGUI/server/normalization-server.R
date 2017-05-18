# read input FCS files, render box #2
observe({
    x <- input$fcsNorm
    if (is.null(x)) return()
    for (i in seq_len(nrow(x)))
        vals$ffsNorm[[i]] <- flowCore::read.FCS(x[[i, "datapath"]])
    output$box2 <- renderUI(box2)
    js$collapse("box_1")
})

# checkboxInput "Normalize to median level of current files"
observe({
    x <- input$box_NormToCurrent
    if (is.null(x) || x == 0) 
        return()
    updateCheckboxInput(session, "box_uploadNormTo", value=FALSE)
    if (length(vals$ffsNorm) == 1) {
        vals$ffsNormTo <- vals$ffsNorm
    } else {
        vals$ffsNormTo <- list(CATALYST::concatFCS(vals$ffsNorm))
    }
    js$collapse("box_2")
})

# checkboxInput "Upload FCS file(s) of beads to normalize to"
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
    } else {
        vals$ffsNormTo <- list(flowCore::read.FCS(input$fcsNormTo$datapath))
    }
    js$collapse("box_2")
})

# use DVS or Beta beads
observe({
    x <- input$select_beads
    if (is.null(x)) return()
    if (x == "dvs" || x == "beta")
        vals$beads <- input$select_beads
})

# selectInput and actionButton for custom beads
output$select_customBeads <- renderUI({
    if (input$select_beads == "custom") 
        selectInput_customBeads(vals$ffsNormTo[[1]])
})

observeEvent(input$button_customBeads, {
    x <- input$input_customBeads
    if (is.null(x)) return()
    chs <- flowCore::colnames(vals$ffsNormTo)
    ms <- gsub("[[:alpha:][:punct:]]", "", chs)
    beadMs <- gsub("[[:alpha:][:punct:]]", "", unlist(x))
    if (sum(beadMs %in% ms) != length(beadMs)) {
        showNotification(ui="Not all bead channels are valid.", 
                         type="error", 
                         closeButton=FALSE)
        return()
    }
    vals$beads <- as.numeric(beadMs)
})

# initialize logical vector of bead indices for each sample
# and list for normalized flowFrame(s), render box #3
observe({
   if (any(vapply(vals$ffsNormTo, is.null, logical(1))))
       return()
    vals$beadInds  <- vector(mode="list", length=length(ffs()))
    vals$beadGates <- vector(mode="list", length=length(ffs()))
    vals$ffsNormed <- vector(mode="list", length=length(vals$ffsNorm))
    output$box3 <- renderUI(box3)
})

# when beads have been selected,
# render gating panel and box #4
observe({
    if (is.null(vals$beads)) return()
    output$box4 <- renderUI(box4)
    output$box_beadGating <- renderUI(
        box_beadGating(samples=smpls(), selected=1))
    js$collapse("box_3")
})

# keep track of currently selected sample
observe({
    x <- input$select_sample
    if (is.null(x)) return()
    vals$smpl <- which(smpls() == x)
})

# list of all input flowFrames
ffs <- reactive({
    if (input$box_NormToCurrent == 1) {
        vals$ffsNorm
    } else if (input$box_uploadNormTo == 1) {
        c(vals$ffsNorm, vals$ffsNormTo)
    } 
})

# next / previous sample buttons
observe({
    toggleState(id="prevSmpl", condition=vals$smpl != 1)
    toggleState(id="nextSmpl", condition=vals$smpl != length(vals$beadInds))
})

observeEvent(input$prevSmpl, { 
    updateSelectInput(session, "select_sample", selected=smpls()[vals$smpl-1]) 
})

observeEvent(input$nextSmpl, { 
    updateSelectInput(session, "select_sample", selected=smpls()[vals$smpl+1]) 
})

# enable gating button only if all gates have been drawn
observe({
    vals$beadGates[[vals$smpl]] <- list(
        input$beadGate1, input$beadGate2, input$beadGate3, input$beadGate4, input$beadGate5)
    test <- !(any(vapply(vals$beadGates[[vals$smpl]], is.null, logical(1))))
    toggleState(id="gateBeads", condition=test)
})

# bead gating
observeEvent(input$gateBeads, {
    selected <- isolate(vals$smpl)
    es <- flowCore::exprs(ffs()[[selected]])
    beadChs <- colnames(es)[isolate(vals$beadCols)]
    dnaCh <- colnames(es)[isolate(vals$dnaCols)][1]
    cf <- isolate(input$input_cfGating)
    
    # get indices of events falling in intersection of all gates
    vals$beadInds[[selected]] <- 
        brushedPoints(brush=vals$beadGates[[selected]][[1]], 
                      allRows=TRUE, df=data.frame(asinh(es/cf)), 
                      xvar=beadChs[1], yvar=dnaCh)[, ncol(es)+1] &
        brushedPoints(brush=vals$beadGates[[selected]][[2]], 
                      allRows=TRUE, df=data.frame(asinh(es/cf)), 
                      xvar=beadChs[2], yvar=dnaCh)[, ncol(es)+1] &
        brushedPoints(brush=vals$beadGates[[selected]][[3]], 
                      allRows=TRUE, df=data.frame(asinh(es/cf)), 
                      xvar=beadChs[3], yvar=dnaCh)[, ncol(es)+1] &
        brushedPoints(brush=vals$beadGates[[selected]][[4]], 
                      allRows=TRUE, df=data.frame(asinh(es/cf)), 
                      xvar=beadChs[4], yvar=dnaCh)[, ncol(es)+1] &
        brushedPoints(brush=vals$beadGates[[selected]][[5]], 
                      allRows=TRUE, df=data.frame(asinh(es/cf)), 
                      xvar=beadChs[5], yvar=dnaCh)[, ncol(es)+1]
})

# valueBox: x/n samples gated
observe({
    n <- length(vals$beadInds)
    x <- n - sum(vapply(vals$beadInds, is.null, logical(1)))
    output$howManyGated <- renderText(paste(x, "/", n))
})

# valueBox: gating yield
observe({
    y <- vals$beadInds[[vals$smpl]]
    if (is.null(y)) { 
        yield <- "0.00"
    } else {
        x <- ffs()[[vals$smpl]]
        yield <- sprintf("%2.2f", sum(y)/nrow(x)*100)
    }
    output$gatingYield <- renderText(yield)
})

# normalization
observe({
    # check that all samples have been gated
    if (any(vapply(vals$beadInds, is.null, logical(1))))
        return()
    n1 <- length(vals$ffsNorm)
    n2 <- length(vals$ffsNormTo)
    if (input$box_NormToCurrent == 1) {
        inds <- which(unlist(vals$beadInds[1:n1]))
    } else {
        inds <- which(unlist(vals$beadInds[(n1+1):(n1+n2)]))
    }

    # get baseline values from ffsNormTo and normalize ffsNorm
    if (n2 > 1) {
        ff <- CATALYST::concatFCS(vals$ffsNormTo)
    } else {
        ff <- vals$ffsNormTo[[1]]
    }
    bl <- getBaseline(ff, vals$beads, inds)
    for (i in seq_along(vals$ffsNorm)) 
        vals$ffsNormed[[i]] <- normCytof(
            vals$ffsNorm[[i]], vals$beads, vals$beadInds[[i]], bl)
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
    
    colnames(smoothedBeads) <- colnames(smoothedBeadsNormed) <- 
        colnames(es)[c(vals$timeCol, vals$beadCols)]
    
    output$plot_smoothedBeads <- renderPlot(
        outPlots2(smoothedBeads, smoothedBeadsNormed, NULL))
})

# ------------------------------------------------------------------------------
# bead vs. dna scatters
# ------------------------------------------------------------------------------
observe({
    beadChs <- vals$beadCols
    dnaCh   <- vals$dnaCols[1]
    if (is.null(beadChs) || is.null(dnaCh)) 
        return()
    (input$button_viewGating)
    es <- flowCore::exprs(ffs()[[vals$smpl]])
    cf <- isolate(input$input_cfGating)
    n  <- isolate(input$input_nGating)

    output$plot_beadScatter1 <- renderPlot(
        plotScatter(es=es, x=beadChs[1], y=dnaCh, cf=cf, n=n))
    output$plot_beadScatter2 <- renderPlot(
        plotScatter(es=es, x=beadChs[2], y=dnaCh, cf=cf, n=n))
    output$plot_beadScatter3 <- renderPlot(
        plotScatter(es=es, x=beadChs[3], y=dnaCh, cf=cf, n=n))
    output$plot_beadScatter4 <- renderPlot(
        plotScatter(es=es, x=beadChs[4], y=dnaCh, cf=cf, n=n))
    output$plot_beadScatter5 <- renderPlot(
        plotScatter(es=es, x=beadChs[5], y=dnaCh, cf=cf, n=n))
})

# ------------------------------------------------------------------------------
# actionButton "Download results"
# ------------------------------------------------------------------------------

output$dwnld_normResults <- downloadHandler(
    filename=function()  
        paste0(Sys.Date(), "_normalization.zip"),
    content=function(file) { 
        dir <- tempdir()
        setwd(dir)
        
        nFiles <- length(vals$ffsNorm)
        nBeads <- numeric(do.call(sum, vals$beadInds))
        isBead <- which(unlist(vals$beadInds))
        
        # generate single plot of concatenated 
        # smoothed beads vs. smoothed normalized beads
        if (nFiles > 1) {
            ff  <- CATALYST::concatFCS(vals$ffsNorm)
            ffN <- CATALYST::concatFCS(vals$ffsNormed)
        } else {
            ff  <- vals$ffsNorm[[1]]
            ffN <- vals$ffsNormed[[1]]
        }
        es  <- flowCore::exprs(ff) [isBead, ]
        esN <- flowCore::exprs(ffN)[isBead, ]
        
        smoothedBeads <- data.frame(
            es[, vals$timeCol],
            vapply(vals$beadCols, function(i)
                stats::runmed(es[, i], 501, "constant"), 
                nBeads))
        
        smoothedBeadsNormed <- data.frame(
            es[, vals$timeCol],
            vapply(vals$beadCols, function(i)
                stats::runmed(esN[, i], 501, "constant"), 
                nBeads))
        
        colnames(smoothedBeads) <- colnames(smoothedBeadsNormed) <- 
            colnames(es)[c(vals$timeCol, vals$beadCols)]
        
        outPlots2(smoothedBeads, smoothedBeadsNormed, out_path=dir)
        
        # write FCS files of normalized data
        # and beads, if "Remove beads" is checked
        first <- TRUE
        for (i in seq_len(nFiles)) {
            if (input$box_removeBeads == 1) {
                if (first) {
                    outNms <- c(rbind(c(
                        gsub(".fcs", "_normed.fcs", input$fcsNorm$name),
                        gsub(".fcs", "_beads.fcs",  input$fcsNorm$name))))
                    first <- FALSE
                }
                suppressWarnings(
                    flowCore::write.FCS(
                        vals$ffsNormed[[i]][!vals$beadInds[[i]], ], 
                        outNms[2*i-1]))
                suppressWarnings(
                    flowCore::write.FCS(
                        vals$ffsNormed[[i]][vals$beadInds[[i]], ],
                        outNms[2*i]))
            } else {
                if (first) {
                    outNms <- c(gsub(".fcs", "_normed.fcs", input$fcsNorm$name))
                    first <- FALSE
                }
                suppressWarnings(
                    flowCore::write.FCS(vals$ffsNormed[[i]], outNms[i]))
            }
        }
        zip(zipfile=file, 
            files=c(outNms, "beads_before_vs_after.pdf")) 
    },
    contentType="application/zip"
)