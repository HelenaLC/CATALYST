# read input FCS, render box #2
observeEvent(input$fcsNorm, {
    # initialize reactive values:
    # list of flowFrame(s) to be normalized,
    # and list of normalized flowFrame(s)
    n <- nrow(input$fcsNorm)
    vals$ffsNorm <- vector("list", n)
    for (i in seq_len(n))
        vals$ffsNorm[[i]] <- flowCore::read.FCS(
            filename=input$fcsNorm[[i, "datapath"]],
            transformation=FALSE,
            truncate_max_range=FALSE)
    vals$ffsNormed <- vector("list", n)
    output$box2 <- renderUI(box2)
    js$collapse("box_1")
})

# checkboxInput "Normalize to median level of current files"
observe({
    x <- input$box_NormToCurrent
    if (is.null(x) || x == 0) 
        return()
    updateCheckboxInput(session, "box_uploadNormTo", value=FALSE)
    # initialize reactive value:
    # list of flowFrame(s) to use as baseline
    if (length(vals$ffsNorm) == 1) {
        vals$ffsNormTo <- vals$ffsNorm[[1]]
    } else {
        vals$ffsNormTo <- CATALYST::concatFCS(vals$ffsNorm)
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
    # initialize reactive value:
    # flowFrame to use as baseline
    vals$ffsNormTo <- flowCore::read.FCS(
        filename=input$fcsNormTo$datapath,
        transformation=FALSE,
        truncate_max_range=FALSE)
    js$collapse("box_2")
})

############################## is this necessary???
smplNmsNorm <- reactive({
    if (input$box_NormToCurrent == 1) {
        input$fcsNorm$name
    } else {
        c(input$fcsNorm$name, input$fcsNormTo$name)
    }
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
        selectInput_customBeads(vals$ffsNormTo)
})

observeEvent(input$button_customBeads, {
    x <- input$input_customBeads
    if (is.null(x)) return()
    # check validity of input bead masses
    chs <- flowCore::colnames(vals$ffsNormTo)
    ms <- gsub("[[:alpha:][:punct:]]", "", chs)
    beadMs <- gsub("[[:alpha:][:punct:]]", "", unlist(x))
    valid <- beadMs %in% ms
    if (sum(valid) != length(beadMs)) {
        showNotification(type="error", closeButton=FALSE,
            ui=paste0("Bead masses ", beadMs[!valid], " are invalid.")) 
        return()
    }
    vals$beads <- as.numeric(beadMs)
})

# get bead, DNA and time column indices
observe({
    if (is.null(vals$beads)) return()
    chs <- flowCore::colnames(ffs()[[vals$smpl]])
    vals$beadCols <- get_bead_cols(chs, isolate(vals$beads))
    vals$dnaCol <- grep("Ir191|Ir193", chs, ignore.case=TRUE)[1]
    vals$timeCol <- grep("time", chs, ignore.case=TRUE)
})

# initialize reactive values:
# logical vector of bead indices for each sample
# and gates, render box #3
observe({
   if (is.null(vals$ffsNormTo)) return()
    vals$beadInds  <- vector("list", length(ffs()))
    vals$beadGates <- vector("list", length(ffs()))
    output$box3 <- renderUI(box3)
})

# when beads have been selected,
# render gating panel and box #4
observe({
    x <- vals$beadCols
    if (is.null(x)) return()
    n <- length(x)
    w <- paste0(100/n, "%")
    output$box4 <- renderUI(box4)
    output$box_beadGating <- renderUI(
        box_beadGating(samples=smplNmsNorm(), selected=1))
    output$beadScatters <- renderUI({
        plotList <- lapply(seq_len(n), function(i) {
            plotOutput(
                outputId=paste0("beadScatter", i), 
                width=w,
                brush=brushOpts(
                    id=paste0("beadGate", i), 
                    resetOnNew=TRUE))
        })
        do.call(tagList, plotList)
    })
    js$collapse("box_3")
})

# keep track of currently selected sample
observe({
    x <- input$select_sample
    if (is.null(x)) return()
    vals$smpl <- which(smplNmsNorm() == x)
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
    updateSelectInput(session, "select_sample", 
        selected=smplNmsNorm()[vals$smpl-1]) 
})

observeEvent(input$nextSmpl, { 
    updateSelectInput(session, "select_sample", 
        selected=smplNmsNorm()[vals$smpl+1]) 
})

# enable gating button only if all gates have been drawn
observe({
    x <- vals$beads
    if (is.null(x)) return()
    selected <- isolate(vals$smpl)
    n <- length(vals$beadCols)
    vals$beadGates[[selected]] <- lapply(seq_len(n), 
        function(i) input[[paste0("beadGate", i)]])
    test <- !(any(vapply(vals$beadGates[[selected]], is.null, logical(1))))
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
    gated <- lapply(seq_len(length(beadChs)), function(i)
        brushedPoints(brush=vals$beadGates[[selected]][[i]], 
            allRows=TRUE, df=data.frame(asinh(es/cf)), 
            xvar=beadChs[i], yvar=dnaCh)[, ncol(es)+1])
    gated <- do.call(cbind, gated)
    vals$beadInds[[selected]] <- apply(gated, 1, function(i)
        sum(i) == length(i))
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

# ------------------------------------------------------------------------------
# normalization
# ------------------------------------------------------------------------------
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
    if (is.null(vals$beadCols) || is.null(vals$dnaCol)) 
        return()
    (input$button_viewGating)
    es <- flowCore::exprs(ffs()[[vals$smpl]])
    cf <- isolate(input$input_cfGating)
    n  <- isolate(input$input_nGating)
    nBeads <- length(vals$beadCols)
    for (i in seq_len(nBeads)) {
        plotId <- paste0("beadScatter", i)
        output[[plotId]] <- renderPlot(
            plotScatter(es=es, x=vals$beadCols[i], y=vals$dnaCol, cf=cf, n=n))
    }
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