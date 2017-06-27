# ==============================================================================
# COMPENSATION
# ==============================================================================

# read input FCS files
ffsComp <- reactive({
    x <- input$fcsComp
    if (is.null(x)) return()
    lapply(seq_len(nrow(x)), function(i) 
        read.FCS(
            filename=x[[i, "datapath"]],
            transformation=FALSE,
            truncate_max_range=FALSE))
})

# expand sidebar once data has been uploaded
output$compensationSidebar1 <- renderUI({
    if (is.null(input$fcsComp)) 
        return()
    compensationSidebar1
})

# ------------------------------------------------------------------------------
# toggle checkboxes
# ------------------------------------------------------------------------------
observe({
    x <- input$box_estSm
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_upldSm", value=FALSE)
})
observe({
    x <- input$box_upldSm
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_estSm", value=FALSE)
})
observe({
    x <- input$box_IDsAsNms
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_upldNms", value=FALSE)
})
observe({
    x <- input$box_upldNms
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_IDsAsNms", value=FALSE)
})

# ------------------------------------------------------------------------------

# render fileInput if checkboxInput == 1
output$inputSm <- renderUI({
    if (input$box_upldSm)
        fileInput("inputSm", NULL, accept=".csv")
})

# get spillover matrix
spillMat <- reactive({
    test <- input$box_upldSm == 1 && !is.null(input$inputSm)
    if (test) {
        read.csv(input$inputSm$datapath, check.names=FALSE, row.names=1)
    } else if (input$box_estSm == 1 && !is.null(dbFrame())) {
        computeSpillmat(x=dbFrame())
    }
})

# make copy for spill adjustment
observe({
    if (is.null(spillMat()))
        return()
    vals$sm <- spillMat()
})

# plotSpillmat()
output$plotSpillmat <- renderPlot({
    if (!is.null(vals$sm)) {
        if (input$box_upldSm) {
            ms <- gsub("[[:punct:][:alpha:]]", "", colnames(vals$sm))
            plotSpillmat(ms, vals$sm)
        } else if (input$box_estSm) {
            plotSpillmat(debaKey(), vals$sm)
        }
    }
})

# compensate input flowFrame(s)
ffsComped <- reactive({
    if (!is.null(vals$sm))
        lapply(ffsComp(), compCytof, vals$sm)
})

# render download button when data has been compensated
output$compensationSidebar2 <- renderUI({
    if (is.null(ffsComped()))
        return()
    compensationSidebar2
})

# ------------------------------------------------------------------------------
# before vs. after compensation scatters
# ------------------------------------------------------------------------------
output$panel_scatters <- renderUI({
    if (is.null(spillMat()))
        return()
    panel_scatters(samples=input$fcsComp$name)
})

# keep track of currently selected sample
selectedSmplComp <- reactive({
    selected <- input$select_smplComp
    if (!is.null(selected))
        match(selected, input$fcsComp$name)
})

observe({
    x <- selectedSmplComp()
    if (is.null(x))
        return()
    choices <- colnames(ffsComp()[[x]])
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
    if (is.null(selectedSmplComp())) return()
    chs <- flowCore::colnames(isolate(vals$ffsComped)[[selectedSmplComp()]])
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
    updateSelectInput(session,
                      inputId="scatterCh1",
                      selected=input$scatterCh2)
    updateSelectInput(session,
                      inputId="scatterCh2",
                      selected=ch1)
})

# textOutput: Spillover of current interaction
output$text_spill <- renderText({
    spill <- vals$sm[ch1(), ch2()]
    if (is.null(spill) || is.na(spill))
        return()
    paste0(sprintf("%.3f", 100*spill), "%")
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

ch1 <- reactive({input$scatterCh1})
ch2 <- reactive({input$scatterCh2})
cfComp <- reactive({input$cfComp})

# default to cofactor 5 if input is invalid
observe({
    if (!is.numeric(input$cfComp) || input$cfComp == 0)
        updateNumericInput(session, inputId="cfComp", value=5)
})

# left scatter (uncompensated)
output$compScatter1 <- renderPlot({ 
    data <- ffsComp()[[selectedSmplComp()]]
    x <- ch1(); y <- ch2()
    if (is.null(data) || x == "" || y == "") return()
    plotScatter(es=exprs(data), n=10e3, x=x, y=y, cf=cfComp())
})
output$text_info1 <- renderText({ 
    data <- ffsComp()[[selectedSmplComp()]]
    x <- ch1(); y <- ch2()
    if (is.null(data) || x == "" || y == "") return()
    text_info(data, isolate(input$cfComp), input$rect1, x, y)
})   

# left scatter (compensated)
output$compScatter2 <- renderPlot({
    data <- ffsComped()[[selectedSmplComp()]]
    x <- ch1(); y <- ch2()
    if (is.null(data) || x == "" || y == "") return()
    plotScatter(es=exprs(data), n=10e3, x=x, y=y, cf=cfComp())
})
output$text_info2 <- renderText({ 
    data <- ffsComped()[[selectedSmplComp()]]
    x <- ch1(); y <- ch2()
    if (is.null(data) || x == "" || y == "") return()
    text_info(data, isolate(input$cfComp), input$rect2, x, y)
})
