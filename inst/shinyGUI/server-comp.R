# ==============================================================================
# COMPENSATION
# ==============================================================================

# read input FCS files
ffsComp <- reactive({
    inFile <- input$fcsComp
    req(inFile)
    lapply(seq_len(nrow(inFile)), function(i) 
        flowCore::read.FCS(
            filename=inFile[[i, "datapath"]],
            transformation=FALSE,
            truncate_max_range=FALSE))
})

# expand sidebar once data has been uploaded
output$compensationSidebar1 <- renderUI({
    req(input$fcsComp)
    compensationSidebar1
})

# toggle checkboxes
observe({
    req(input$box_estSm == 1)
    updateCheckboxInput(session, "box_upldSm", value=FALSE)
})
observe({
    req(input$box_upldSm == 1)
    updateCheckboxInput(session, "box_estSm", value=FALSE)
})

# ------------------------------------------------------------------------------
# spillover matrix selection & plotting
# ------------------------------------------------------------------------------

# render fileInput if checkboxInput == 1
output$inputSm <- renderUI({
    req(input$box_upldSm == 1)
    fileInput("inputSm", NULL, accept=".csv")
})

# get spillover matrix
spillMat <- reactive({
    req(ffsComp())
    x <- input$box_upldSm
    if (!is.null(x) && x == 1 && !is.null(input$inputSm)) {
        read.csv(input$inputSm$datapath, check.names=FALSE, row.names=1)
    } else {
        x <- input$box_estSm
        if (!is.null(x) && x == 1 && !is.null(dbFrame()))
            computeSpillmat(x=dbFrame())
    }
})

# store original spillover matrix
observeEvent(spillMat(), {
    vals$sm <- spillMat()
})

# plotSpillmat()
output$plotSpillmat <- renderPlot({
    req(vals$sm)
    if (input$box_upldSm) {
        ms <- gsub("[[:punct:][:alpha:]]", "", colnames(vals$sm))
        CATALYST::plotSpillmat(bc_ms=ms, SM=vals$sm)
    } else if (input$box_estSm) {
        CATALYST::plotSpillmat(bc_ms=debaKey(), SM=vals$sm)
    }
})

# compensate input flowFrame(s)
ffsComped <- reactive({
    req(vals$sm)
    lapply(ffsComp(), CATALYST::compCytof, vals$sm)
})

# render download button when data has been compensated
output$compensationSidebar2 <- renderUI({
    req(ffsComped())
    compensationSidebar2
})

# ------------------------------------------------------------------------------
# before vs. after compensation scatters
# ------------------------------------------------------------------------------
output$panel_scatters <- renderUI({
    req(spillMat())
    panel_scatters(samples=input$fcsComp$name)
})

# keep track of currently selected sample
selectedSmplComp <- reactive({
    req(input$select_smplComp)
    match(input$select_smplComp, input$fcsComp$name)
})

observe({
    selected <- selectedSmplComp()
    req(selected)
    choices <- colnames(ffsComp()[[selected]])
    ms <- gsub("[[:alpha:][:punct:]]", "", choices)
    ch1 <- choices[choices == input$scatterCh1]
    ch2 <- choices[choices == input$scatterCh2]
    if (length(ch1) != 1) {
        m <- gsub("[[:alpha:][:punct:]]", "", ch1)
        ch1 <- choices[match(m, ms)]
    }
    if (length(ch2) != 1) {
        m <- gsub("[[:alpha:][:punct:]]", "", ch1)
        ch1 <- choices[match(m, ms)]
    }
    updateSelectInput(session,
        inputId="scatterCh1",
        choices=choices,
        selected=ch1)
    updateSelectInput(session,
        inputId="scatterCh2",
        choices=choices,
        selected=ch2)
})

# ------------------------------------------------------------------------------
# next / previous sample buttons
# ------------------------------------------------------------------------------
observe({
    n <- length(ffsComp())
    toggleState(id="prev_smplComp", condition=selectedSmplComp() != 1)
    toggleState(id="next_smplComp", condition=selectedSmplComp() != n)
})
observeEvent(input$prev_smplComp, { 
    updateSelectInput(session, 
        inpudId="select_smplComp", 
        selected=input$fcsComp$name[selectedSmplComp()-1]) 
})
observeEvent(input$next_smplComp, { 
    updateSelectInput(session, 
        inputId="select_smplComp", 
        selected=input$fcsComp$name[selectedSmplComp()+1]) 
})

# update channel selection choices upon sample change
observe({
    selected <- selectedSmplComp()
    req(selected)
    chs <- flowCore::colnames(isolate(vals$ffsComped)[[selected]])
    # default to currently selected masses 
    ms <- gsub("[[:punct:][:alpha:]]", "", chs)
    currentMs <- sapply(
        c(isolate(input$scatterCh1), isolate(input$scatterCh2)), 
        function(i) gsub("[[:punct:][:alpha:]]", "", i))
    inds <- match(currentMs, ms)
    # if non-existent, plot first 2 mass channels
    if (any(is.na(inds))) 
        inds <- which(!is.na(as.numeric(ms)))[1:2]
    updateSelectInput(session, 
        inputId="scatterCh1", 
        choices=chs, 
        selected=chs[inds[1]])
    updateSelectInput(session, 
        inputId="scatterCh2", 
        choices=chs, 
        selected=chs[inds[2]])
})

# actionButton: Swap x- and y-axis
observeEvent(input$flipAxes, {
    ch1 <- ch1()
    updateSelectInput(session, inputId="scatterCh1", selected=input$scatterCh2)
    updateSelectInput(session, inputId="scatterCh2", selected=ch1)
})

# textOutput: Spillover of current interaction
output$text_spill <- renderText({
    ch1 <- ch1()
    ch2 <- ch2()
    if (ch1() %in% rownames(vals$sm) && ch2() %in% colnames(vals$sm)) {
        spill <- vals$sm[ch1, ch2]
        paste0(sprintf("%.3f", 100*spill), "%")
    }
})

# actionButton: Adjust spill of current interaction
observe(toggleState(id="adjustSpill", condition=is.numeric(input$newSpill)))
observeEvent(input$adjustSpill, {
    vals$sm[ch1(), ch2()] <- input$newSpill/100
})

# actionButtons revert current / all spill adjustments
observeEvent(input$revert, {
    vals$sm[ch1(), ch2()] <- spillMat()[ch1(), ch2()]
})
observeEvent(input$revertAll, {
    vals$sm <- spillMat()
})

# ------------------------------------------------------------------------------
# before vs. after compensation scatters
# ------------------------------------------------------------------------------

ch1 <- reactive(input$scatterCh1)
ch2 <- reactive(input$scatterCh2)
cfComp <- reactive(input$cfComp)

# default to cofactor 5 if input is invalid
observe({
    req(!is.numeric(input$cfComp), input$cfComp == 0)
    updateNumericInput(session, inputId="cfComp", value=5)
})

# left scatter (uncompensated)
output$compScatter1 <- renderPlot({ 
    data <- ffsComp()[[selectedSmplComp()]]
    x <- ch1(); y <- ch2()
    req(!is.null(data), x != "", y != "")
    CATALYST:::plotScatter(es=flowCore::exprs(data), x=x, y=y, cf=cfComp())
})
output$text_info1 <- renderText({ 
    data <- ffsComp()[[selectedSmplComp()]]
    x <- ch1(); y <- ch2()
    req(!is.null(data), x != "", y != "")
    text_info(data, isolate(input$cfComp), input$rect1, x, y)
})   

# left scatter (compensated)
output$compScatter2 <- renderPlot({
    data <- ffsComped()[[selectedSmplComp()]]
    x <- ch1(); y <- ch2()
    req(!is.null(data), x != "", y != "")
    CATALYST:::plotScatter(es=flowCore::exprs(data), x=x, y=y, cf=cfComp())
})
output$text_info2 <- renderText({ 
    data <- ffsComped()[[selectedSmplComp()]]
    x <- ch1(); y <- ch2()
    req(!is.null(data), x != "", y != "")
    text_info(data, isolate(input$cfComp), input$rect2, x, y)
})

# ------------------------------------------------------------------------------
# download handlers
# ------------------------------------------------------------------------------

output$dwnld_debaPlots <- downloadHandler(
    filename=function() { "yield_event_plots.zip" },
    content =function(file) { 
        dbFrame <- dbFrame()
        tmpdir <- tempdir()
        setwd(tmpdir)
        CATALYST::plotYields(
            x=dbFrame, 
            which=yieldPlotChoices(), 
            out_path=tmpdir)
        CATALYST::plotEvents(
            x=dbFrame, 
            out_path=tmpdir, 
            n_events=250)
        zip(zipfile=file,
            files=paste0(c("yield_plot", "event_plot"), ".pdf")) },
    contentType="application/zip")

output$dwnld_comped <- downloadHandler(
    filename=function() { 
        paste0(Sys.Date(), "_compensation.zip")
    },
    content =function(file) { 
        tmpdir <- tempdir()
        setwd(tmpdir)
        outNms <- file.path(tmpdir, gsub(".fcs", "_comped.fcs", 
            input$fcsComp, ignore.case=TRUE))
        ffs <- ffsComped()
        if (input$box_setToZero) {
            lapply(seq_along(ffs), function(i) {
                flowCore::exprs(ffs[[i]])[flowCore::exprs(ffs[[i]]) < 0] <- 0
                suppressWarnings(flowCore::write.FCS(ffs[[i]], outNms[i])) 
            })
        } else {
            lapply(seq_along(ffs), function(i) {
                suppressWarnings(flowCore::write.FCS(ffs[[i]], outNms[i]))
            })
        }
        zip(zipfile=file, files=gsub(tmpdir, "", outNms)) 
    },
    contentType="application/zip")

output$dwnld_spillMat <- downloadHandler(
    filename=function() { 
        "spillMat.csv"
    },
    content=function(file) { 
        write.csv(vals$sm, file) 
    })