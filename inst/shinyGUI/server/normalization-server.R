# ==============================================================================
# NORMALIZATION
# ==============================================================================

# read input FCS
ffsNorm <- reactive({
    x <- input$fcsNorm
    if (is.null(x)) return()
    lapply(seq_len(nrow(x)), function(i)
        flowCore::read.FCS(
            filename=x[[i, "datapath"]],
            transformation=FALSE,
            truncate_max_range=FALSE))
})

output$box2 <- renderUI({
    if (!is.null(ffsNorm())) {
        js$collapse("box1")
        box2
    }
})
output$box3 <- renderUI({
    if (!is.null(ffNormTo())) {
        js$collapse("box2")
        box3
    }
})
output$box4 <- renderUI({
    if (!is.null(beadCols())) {
        js$collapse("box3")
        box4
    }
})

# toggle checkboxes
observe({
    x <- input$box_NormToCurrent
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_uploadNormTo", value=FALSE)
})
observe({
    x <- input$box_uploadNormTo
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_NormToCurrent", value=FALSE)
})

# render fileInput if checkboxInput == 1
output$input_NormTo <- renderUI({
    x <- input$box_uploadNormTo
    if (is.null(x) || x == 0)
        return()
    fileInput(
        inputId="fcsNormTo", 
        label=NULL, 
        accept=".fcs")
})

# get flowFrame to use as baseline
ffNormTo <- reactive({
    x <- input$box_NormToCurrent
    if (!is.null(x) && x == 1) {
        # use input flowFrame(s)
        if (length(ffsNorm()) == 1) {
            ffsNorm()[[1]]
        } else {
            CATALYST::concatFCS(ffsNorm())
        }
    } else {
        x <- input$box_uploadNormTo
        if (!is.null(x) && x == 1 && !is.null(input$fcsNormTo)) {
            # read FCS of beads to normalize to
            read.FCS(
                filename=input$fcsNormTo$datapath,
                transformation=FALSE,
                truncate_max_range=FALSE)
        }
    }
})

# get custom beads from textInput
observeEvent(input$button_customBeads, {
    x <- input$input_customBeads
    if (is.null(x)) return()
    # check validity of input bead masses
    chs <- flowCore::colnames(ffNormTo())
    ms <- gsub("[[:alpha:][:punct:]]", "", chs)
    beadMs <- gsub("[[:alpha:][:punct:]]", "", unlist(x))
    valid <- beadMs %in% ms
    if (sum(valid) != length(beadMs)) {
        showNotification(type="error", closeButton=FALSE,
            ui=paste0("Bead masses ", beadMs[!valid], " are invalid.")) 
        return()
    }
    vals$customBeads <- as.numeric(beadMs)
})

# get beads: "dvs", "beta", or custom
beads <- reactive({
    x <- input$select_beads
    if (is.null(x)) return()
    if (x == "dvs" || x == "beta") {
        input$select_beads
    } else if (!is.null(vals$customBeads)) {
        vals$customBeads
    }
})

# selectInput and actionButton for custom beads
output$select_customBeads <- renderUI({
    if (input$select_beads == "custom") 
        selectInput_customBeads(ffNormTo())
})

# get current channels, bead, DNA and time column indices
chs <- reactive({
    if (!is.null(ffNormTo()))
        return()
    flowCore::colnames(ffNormTo())
})
beadCols <- reactive({
    if (is.null(beads())) return()
    get_bead_cols(chs(), beads())
})
dnaCol <- reactive({
    if (is.null(beads())) return()
    grep("Ir191|Ir193", chs(), ignore.case=TRUE)[1]
})
timeCol <- reactive({
    if (is.null(beads())) return()
    grep("time", chs(), ignore.case=TRUE)
})

# when beads have been selected,
# render gating panel and box #4
# observe({
#     x <- vals$beadCols
#     if (is.null(x)) return()
#     n <- length(x)
#     w <- paste0(100/n, "%")
#     output$box4 <- renderUI(box4)
#     output$box_beadGating <- renderUI(
#         box_beadGating(samples=smplNmsNorm(), selected=1))
#     output$beadScatters <- renderUI({
#         plotList <- lapply(seq_len(n), function(i) {
#             plotOutput(
#                 outputId=paste0("beadScatter", i), 
#                 width=w,
#                 brush=brushOpts(
#                     id=paste0("beadGate", i), 
#                     resetOnNew=TRUE))
#         })
#         do.call(tagList, plotList)
#     })
#     js$collapse("box3")
# })

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

smplNmsNorm <- reactive({
    if (input$box_NormToCurrent) {
        input$fcsNorm$name
    } else if (input$box_upldNormTo) {
        input$fcsNormTo$name
    }
})

# keep track of currently selected sample
selectedSmplNorm <- reactive({
    x <- input$select_sample
    if (is.null(x)) return()
    which(smplNmsNorm() == x)
})

beadGates <- reactive({
    x <- beadCols()
    if (is.null(x)) 
        return()
    n <- length(x)
    vector("list", n)
})
beadInds <- reactive({
    x <- beadCols()
    if (is.null(x)) 
        return()
    n <- length(x)
    vector("list", n)
})

cfNorm <- reactive({
    if (!is.null(input$input_cfGating))
        input$input_cfGating
})

# enable gating button only if all gates have been drawn
# on currently selected sample
observe({
    x <- selectedSmplComp()
    if (is.null(x)) return()
    test <- !(any(vapply(beadGates()[[x]], is.null, logical(1))))
    toggleState(id="gateBeads", condition=test)
})

# get indices of events falling in intersection of all gates
observeEvent(input$gateBeads, {
    x <- selectedSmplNorm()
    if (is.null(x)) return()
    n <- length(beadCols())
    beadGates()[[x]] <- lapply(seq_len(n), 
        function(i) input[[paste0("beadGate", x)]])
    es <- flowCore::exprs(ffsNormTo()[[x]])
    chs <- colnames(es)
    gated <- lapply(seq_len(n), function(i)
        brushedPoints(
            brush=beadGates()[[x]],
            allRows=TRUE,
            df=data.frame(asinh(es/cfNorm())),
            xvar=chs[beadCols()[i]],
            yvar=chs[dnaCol()])[, ncol(es)+1])
    gated <- do.call(cbind, gated)
    beadInds()[[x]] <- apply(gated, 1, function(i) sum(i) == length(i))
})

# valueBox: x/n samples gated
observe({
    n <- length(beadCols())
    x <- n - sum(vapply(beadInds(), is.null, logical(1)))
    output$howManyGated <- renderText(paste(x, "/", n))
})

# valueBox: gating yield
observe({
    if (is.null(beadInds()))
    y <- beadInds()[[selectedSmplNorm()]]
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
    n1 <- length(ffsNorm())
    n2 <- length(ffsNormTo())
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