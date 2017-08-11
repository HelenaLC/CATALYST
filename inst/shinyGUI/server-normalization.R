# ==============================================================================
# NORMALIZATION
# ==============================================================================

# read input FCS
ffsNorm <- reactive({
    if (vals$keepDataConcat) {
        list(ffConcat())
    } else {
        req(input$fcsNorm)
        n <- nrow(input$fcsNorm)
        # check validity of input FCS files
        valid <- check_FCS_fileInput(input$fcsNorm, n)
        if (!valid) return()
        lapply(seq_len(n), function(i)
            flowCore::read.FCS(
                filename=input$fcsNorm[[i, "datapath"]],
                transformation=FALSE,
                truncate_max_range=FALSE))
    }
})

output$box2 <- renderUI({
    req(ffsNorm())
    js$collapse("box_1")
    box2
})
output$box3 <- renderUI({
    req(input$box_normToCurrent)
    js$collapse("box_2")
    box3
})
output$box4 <- renderUI({
    req(beadCols())
    js$collapse("box_3")
    box4
})

# ------------------------------------------------------------------------------
# Normalization baseline
# ------------------------------------------------------------------------------

# toggle checkboxes
observe({
    req(input$box_normToCurrent == 1)
    updateCheckboxInput(session, "box_uploadNormTo", value=FALSE)
})
observe({
    req(input$box_uploadNormTo == 1)
    updateCheckboxInput(session, "box_normToCurrent", value=FALSE)
})

# render fileInput for FCS of reference beads
output$input_NormTo <- renderUI({
    req(input$box_uploadNormTo == 1)
    fileInput(
        inputId="fcsNormTo", 
        label=NULL, 
        accept=".fcs")
})

# get normalization baseline
baseline <- reactive({
    req(beadCols())
    x <- input$box_uploadNormTo
    if (!is.null(x) && x && !is.null(input$fcsNormTo)) {
        ref_ff <- flowCore::read.FCS(
            filename=input$fcsNormTo$datapath,
            transformation=FALSE,
            truncate_max_range=FALSE)
        colMeans(flowCore::exprs(ref_ff)[, beadCols()])
    } else {
        x <- input$box_normToCurrent
        if (!is.null(x) && x &&
                !(any(vapply(vals$beadInds, is.null, logical(1))))) 
            colMeans(matrix(do.call(rbind, lapply(seq_along(ffsNorm()), 
                function(i) colMeans(flowCore::exprs(
                    ffsNorm()[[i]][vals$beadInds[[i]], beadCols()])))), 
                ncol=length(beadCols())))
    }
})

# ------------------------------------------------------------------------------
# Bead selection
# ------------------------------------------------------------------------------
# selectInput and actionButton for custom beads
output$select_customBeads <- renderUI({
    if (input$select_beads == "custom") 
        selectInput_customBeads(ffsNorm()[[1]])
})

# get custom beads from textInput
observeEvent(input$button_customBeads, {
    req(input$input_customBeads)
    # check validity of input bead masses
    ms <- gsub("[[:alpha:][:punct:]]", "", chs())
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
    if (!is.null(x) && (x == "dvs" || x == "beta")) {
        x
    } else if (!is.null(vals$customBeads)) {
        vals$customBeads
    }
})

# get current channels, bead, DNA, time and length column indices
chs <- reactive({
    req(ffsNorm())
    colnames(ffsNorm()[[1]])
})
beadCols <- reactive({
    req(beads())
    CATALYST:::get_bead_cols(chs(), beads())
})
dnaCol <- reactive({
    req(chs())
    grep("Ir191|Ir193", chs(), ignore.case=TRUE)[1]
})
timeCol <- reactive({
    req(chs())
    grep("time", chs(), ignore.case=TRUE)
})
lengthCol <- reactive({
    req(chs())
    grep("length", chs(), ignore.case=TRUE)
})

# next / previous sample buttons
observe({
    n <- length(ffsNorm())
    toggleState(id="prevSmplGating", condition=selectedSmplGating() != 1)
    toggleState(id="nextSmplGating", condition=selectedSmplGating() != n)
})
observeEvent(input$prevSmplGating, { 
    updateSelectInput(session, "selectSmplGating", 
        selected=smplNmsNorm()[selectedSmplGating()-1]) 
})
observeEvent(input$nextSmplGating, { 
    updateSelectInput(session, "selectSmplGating", 
        selected=smplNmsNorm()[selectedSmplGating()+1]) 
})
observe({
    n <- length(ffsNorm())
    toggleState(id="prevSmplMhl", condition=selectedSmplMhl() != 1)
    toggleState(id="nextSmplMhl", condition=selectedSmplMhl() != n)
})
observeEvent(input$prevSmplMhl, { 
    updateSelectInput(session, "selectSmplMhl", 
        selected=smplNmsNorm()[selectedSmplMhl()-1]) 
})
observeEvent(input$nextSmplMhl, { 
    updateSelectInput(session, "selectSmplMhl", 
        selected=smplNmsNorm()[selectedSmplMhl()+1]) 
})

# keep track of currently selected sample
smplNmsNorm <- reactive({
    if (isTRUE(vals$keepDataConcat)) {
        smplNmConcat()
    } else {
        req(input$fcsNorm)
        input$fcsNorm$name
    }
})
selectedSmplGating <- reactive({
    x <- input$selectSmplGating
    if (is.null(x)) return(1)
    which(smplNmsNorm() == x)
})
selectedSmplMhl <- reactive({
    x <- input$selectSmplMhl
    if (is.null(x)) return(1)
    which(smplNmsNorm() == x)
})

# ------------------------------------------------------------------------------
# Bead gating
# ------------------------------------------------------------------------------
# render gating panel when beads have been selected
output$box_beadGating <- renderUI({
    req(beadCols(), smplNmsNorm())
    box_beadGating(samples=smplNmsNorm())
})

# initialize plot list
output$beadScatters <- renderUI({
    req(beadCols())
    n <- length(beadCols())
    w <- paste0(100/n-1, "%")
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

observe({
    req(beadCols(), dnaCol())
    for (i in seq_along(beadCols())) {
        local({
            j <- i
            plotId <- paste0("beadScatter", j)
            output[[plotId]] <- renderPlot(
                CATALYST:::plotScatter(
                    es=exprs(ffsNorm()[[selectedSmplGating()]]), 
                    x=beadCols()[j], y=dnaCol(), cf=5))
        })
    }
})
  
# initialize list of gates, bead indices,
# and vector of Mahalanobis distance cutoffs
observeEvent(beadCols(), {
    n <- length(ffsNorm())
    vals$beadGates =vector("list", n)
    vals$beadInds  =vector("list", n)
    vals$mhlCutoffs=rep(0, n)
})

# enable gating button only if all gates have been drawn
observe({
    test <- !any(vapply(seq_along(beadCols()), function(i)
        is.null(input[[paste0("beadGate", i)]]), logical(1)))
    toggleState(id="gateBeads", condition=test)
})

# bsButton: "Gate"
observeEvent(input$gateBeads, {
    x <- selectedSmplGating()
    n <- length(beadCols())
    es <- flowCore::exprs(ffsNorm()[[x]])
    # get indices of events falling in intersection of all gates
    vals$beadGates[[x]] <- lapply(seq_len(n), 
        function(i) input[[paste0("beadGate", i)]])
    gated <- lapply(seq_len(n), function(i)
        brushedPoints(
            brush=vals$beadGates[[x]][[i]],
            allRows=TRUE,
            df=data.frame(asinh(es/5)),
            xvar=chs()[beadCols()[i]],
            yvar=chs()[dnaCol()])[, ncol(es)+1])
    gated <- do.call(cbind, gated)
    vals$beadInds[[x]] <- apply(gated, 1, function(i) sum(i) == length(i))
})

# normalize input flowFrames when all samples have been gated
ffsNormed <- reactive({
    # check that all samples have been gated
    req(vals$beadInds, !any(vapply(vals$beadInds, is.null, logical(1))))
    lapply(seq_along(ffsNorm()), function(i) {
        es <- flowCore::exprs(ffsNorm()[[i]])
        # get bead slopes and linearly interpolate at non-bead events
        bead_es <- es[vals$beadInds[[i]], beadCols()]
        bead_ts <- es[vals$beadInds[[i]], timeCol()]
        bead_slopes <- rowSums(bead_es*baseline())/rowSums(bead_es^2)
        slopes <- approx(bead_ts, bead_slopes, es[, timeCol()], rule=2)$y
        # normalize
        ex <- c(timeCol(), lengthCol())
        ffNormed <- ffsNorm()[[i]]
        flowCore::exprs(ffNormed)[, -ex] <- 
            flowCore::exprs(ffNormed)[, -ex] * slopes
        ffNormed
    })
})

# compute Mahalanobis distances from identified beads
mhlDists <- reactive({
    req(input$box_removeBeads == 1,
        # check that all samples have been gated
        !(any(vapply(vals$beadInds, is.null, logical(1)))))
    lapply(seq_along(ffsNorm()), function(i) {
        es <- flowCore::exprs(ffsNorm()[[i]])
        es_t <- asinh(es[, beadCols()]/5)
        sqrt(stats::mahalanobis(x=es_t, 
            center=colMeans(es_t[vals$beadInds[[i]], ]), 
            cov=stats::cov(es_t[vals$beadInds[[i]], ])))
    })
})

# ------------------------------------------------------------------------------
# valueBox: "x/n samples gated"
# ------------------------------------------------------------------------------
output$howManyGated <- renderText({
    if (is.null(beadCols())) 
        return()
    n <- length(ffsNorm())
    if (is.null(vals$beadInds)) {
        y <- 0
    } else {
        y <- n - sum(vapply(vals$beadInds, is.null, logical(1)))
    }
    paste(y, "/", n)
})

# ------------------------------------------------------------------------------
# # valueBox: "gating yield"
# ------------------------------------------------------------------------------
output$gatingYield <- renderText({
    if (is.null(beadCols())) 
        return()
    selected <- selectedSmplGating()
    inds <- vals$beadInds[[selected]]
    if (is.null(inds)) {
        yield <- "0.00"
    } else {
        ff <- ffsNorm()[[selected]]
        yield <- sprintf("%2.2f", sum(inds)/nrow(ff)*100)
    }
    yield
})

# ------------------------------------------------------------------------------
# Plot smoothed beads vs. time
# ------------------------------------------------------------------------------
# render UI once ffsNorm have been normalized
output$box_smoothedBeads <- renderUI({
    req(ffsNormed())
    js$collapse("beadGating")
    box_smoothedBeads
})

smoothedBeads <- reactive({
    req(ffsNormed())
    selected <- selectedSmplGating()
    beadInds <- vals$beadInds[[selected]]
    if (sum(beadInds) == 0) {
        showNotification("No beads have been gated for this sample. 
            Please consider re-gating.", type="error", closeButton=FALSE)
    }
    es <- flowCore::exprs(ffsNorm()[[selected]])[beadInds, ]
    normedEs <- flowCore::exprs(ffsNormed()[[selected]])[beadInds, ]
    
    n <- sum(beadInds)
    smoothed <- data.frame(
        es[, timeCol()],
        vapply(beadCols(), function(i)
            stats::runmed(es[, i], 501, "constant"), numeric(n)))
    smoothedNormed <- data.frame(
        es[, timeCol()],
        vapply(beadCols(), function(i)
            stats::runmed(normedEs[, i], 501, "constant"), numeric(n)))
    colnames(smoothed) <- colnames(smoothedNormed) <- 
        chs()[c(timeCol(), beadCols())]
    
    p1 <- CATALYST:::plotSmoothed(smoothed, "Smoothed beads")
    p2 <- CATALYST:::plotSmoothed(smoothedNormed, "Smoothed normalized beads")
    CATALYST:::arrangeSmoothed(p1, p2, shiny=TRUE)
})

output$plot_smoothedBeads <- renderPlot(grid.arrange(smoothedBeads()))

# ------------------------------------------------------------------------------
# bead removal
# ------------------------------------------------------------------------------
# render sample selection, sliderInput & actionButton
output$mhlCutoffNormUI <- renderUI({
    req(input$box_removeBeads == 1, mhlDists())
    x <- selectedSmplMhl()
    mhlCutoffNormUI(
        samples=smplNmsNorm(), 
        selected=x,
        maxMhlDist=ceiling(max(mhlDists()[[x]])))
})

# render beads vs. beads panel
output$box_beadRemoval <- renderUI({
    req(input$box_removeBeads == 1, mhlDists())
    js$collapse("smoothedBeads")
    box_beadRemoval
})

# apply cutoff to currently selected sample
observeEvent(input$applyMhlCutoffNorm, {
    x <- selectedSmplMhl()
    vals$mhlCutoffs[x] <- input$mhlCutoffNorm
})

# render beads vs. beads plot color coded according to 
# Mahalanobis distance from identified beads
output$beadsVsBeads <- renderPlot({
    req(mhlDists())
    x <- selectedSmplMhl()
    es <- asinh(exprs(ffsNorm()[[x]])[, beadCols()]/5)
    tmp <- CATALYST:::get_axes(es, 5)
    tcks <- tmp[[1]]
    labs <- tmp[[2]]
    CATALYST:::plotBeadsVsBeads(
        es, mhlDists()[[x]], vals$mhlCutoffs[x], tcks, labs)
})

# ------------------------------------------------------------------------------
# enable "Go to debarcoding" and downloadButton 
# once data has been normalized & (optionally) 
# Mahalanobis cutoffs have been applied
observe({
    req(beadCols())
    x <- input$box_removeBeads
    test <- x == 0 && !is.null(ffsNormed()) ||
        x == 1 && !any(vals$mhlCutoffs == 0)
    toggleState(id="goToComp", condition=test)
    toggleState(id="dwnld_normResults", condition=test)
})

# bsButton "Go to debarcoding": propagate data 
# & hide FCS fileInput from compensation tab
observeEvent(input$goToComp, {    
    vals$keepDataNorm <- TRUE
    shinyjs::hide(id="fcsComp")
    updateTabItems(session, inputId="tabs", selected="compensation")
})

# ------------------------------------------------------------------------------
# download handler
# ------------------------------------------------------------------------------
output$dwnld_normResults <- downloadHandler(
    filename=function()  
        paste0(format(Sys.Date(), "%y%m%d"), "_normalization.zip"),
    content=function(file) { 
        dir <- tempdir()
        setwd(dir)
        # plot smoothed beads before vs. after normalization
        gt <- smoothedBeads()
        ggsave("beads_before_vs_after.pdf", width=12, height=8, plot=gt)
        # write FCS files of normalized data,
        # beads, and removed events
        raw <- ffsNorm()
        normed <- ffsNormed()
        nms <- gsub(".fcs", "", smplNmsNorm(), ignore.case=TRUE)
        outNms <- rbind(
            paste0(nms, "_normed.fcs"),
            paste0(nms, "_beads.fcs"),
            paste0(nms, "_removed.fcs"))
        if (input$box_removeBeads) {
            dists <- mhlDists()
            sapply(seq_along(raw), function(i) {
                beads <- vals$beadInds[[i]]
                removed <- dists[[i]] < vals$mhlCutoffs[i]
                suppressWarnings(flowCore::write.FCS(
                    x=normed[[i]][!removed, ], 
                    filename=outNms[1, i]))
                suppressWarnings(flowCore::write.FCS(
                    x=raw[[i]][beads, ],  
                    filename=outNms[2, i]))
                suppressWarnings(flowCore::write.FCS(
                    x=normed[[i]][removed, ], 
                    filename=outNms[3, i]))
            })
        } else {
            sapply(seq_along(raw), function(i) {
            beads <- vals$beadInds[[i]]
            suppressWarnings(flowCore::write.FCS(
                x=normed[[i]],
                filename=outNms[1, i]))
            suppressWarnings(flowCore::write.FCS(
                x=raw[[i]][beads, ],  
                filename=outNms[2, i]))
            })
        }
        zip(zipfile=file, 
            files=c("beads_before_vs_after.pdf", outNms)) 
    },
    contentType="application/zip"
)