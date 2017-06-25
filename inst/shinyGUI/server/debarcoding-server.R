# ==============================================================================
# DEBARCODING
# ==============================================================================

# read input FCS
ffDeba <- reactive(flowCore::read.FCS(input$fcsDeba$datapath))

# expand sidebar
output$debarcodingSidebar1 <- renderUI({
    if (is.null(input$fcsDeba)) return()
    debarcodingSidebar1
})

# ------------------------------------------------------------------------------
# checkboxInputs 
# "Upload barcoding scheme (CSV)" &
# "Select single-positive channels"
# ------------------------------------------------------------------------------

# toggle checkboxes
observe({
    x <- input$boxUploadCsv
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "boxSelectBcChs", value=FALSE)
})
observe({
    x <- input$boxSelectBcChs
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "boxUploadCsv", value=FALSE)
})

# boxUploadCsv -> render fileInput
output$debaSchemeCsv <- renderUI({
    x <- input$boxUploadCsv
    if (is.null(x) || x == 0) return()
    fileInput(inputId="debaSchemeCsv", label=NULL, accept=".csv")
})

# boxSelectBcChs -> render selectInput
output$selectBcChs <- renderUI({
    x <- input$boxSelectBcChs
    if (is.null(x) || x == 0) return()
    selectInput(
        inputId="input_bcChs", label=NULL, 
        choices=flowCore::colnames(ffDeba()), 
        multiple=TRUE, selectize=FALSE, size=12)
})

# get debarcoding scheme
debaKey <- reactive({
    if (!is.null(input$input_bcChs)) {
        gsub("[[:alpha:][:punct:]]", "", input$input_bcChs)
    } else if (!is.null(input$debaSchemeCsv)) {
        read.csv(input$debaSchemeCsv$datapath, check.names=FALSE, row.names=1)
    } else {
        NULL
    }
})

# ------------------------------------------------------------------------------
# bsButton: "Debarcode"
# ------------------------------------------------------------------------------
observeEvent(input$buttonDebarcode, {
    if (is.null(debaKey())) {
        showNotification("Please upload a barcoding scheme or select 
            single-positive channels.", type="error", closeButton=FALSE)
        return()
    }
    showNotification(h5("Assigning preliminary IDs..."), 
        duration=NULL, closeButton=FALSE, id="msg", type="message")
    vals$dbFrame1 <- assignPrelim(ffDeba(), debaKey())
    removeNotification(id="msg")
    output$debarcodingSidebar2 <- renderUI(debarcodingSidebar2)
})

# ------------------------------------------------------------------------------
# toggle checkboxes
observe({
    x <- input$box_estCutoffs
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_adjustCutoff", value=FALSE)
    updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
})
observe({
    x <- input$box_adjustCutoff
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_estCutoffs",   value=FALSE)
    updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
})
observe({
    x <- input$box_globalCutoff
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_estCutoffs",   value=FALSE)
    updateCheckboxInput(session, "box_adjustCutoff", value=FALSE)
})

# ------------------------------------------------------------------------------
# cutoff estimation
# ------------------------------------------------------------------------------
observeEvent(input$box_estCutoffs, {
    vals$dbFrame2 <- estCutoffs(vals$dbFrame1)
})

# ------------------------------------------------------------------------------
# cutoff adjustment
# ------------------------------------------------------------------------------

# toggle selectInput, numericInput & bsButton
observe({
    x <- input$box_adjustCutoff
    toggle(id="adjustCutoff", condition=x)
})

# selectInput choices
adjustCutoffChoices <- reactive({
    x <- dbFrame()
    if (is.null(x)) {
        NULL
    } else {
        rownames(debaKey())
    }
})

# renderUI if checkboxInput == 1
output$adjustCutoffUI <- renderUI({
        x <- adjustCutoffChoices() == input$select_adjustCutoff
        adjustCutoffUI(
            dbFrame=vals$dbFrame2, 
            choices=adjustCutoffChoices(), 
            selected=x)
})

# synchronize selectInput w/ yield plot
observe({
    x <- input$select_yieldPlot
    if (is.null(x) || x == 0) return()
    updateSelectInput(session, "select_adjustCutoff", selected=x)
})

# adjust cutoff upon bsButton click
observeEvent(input$button_adjustCutoff, {
    x <- adjustCutoffChoices() == input$select_adjustCutoff
    sep_cutoffs(vals$dbFrame2)[x] <- input$input_adjustCutoff
})

# ------------------------------------------------------------------------------
# global separation cutoff
# ------------------------------------------------------------------------------

# toggle numericInput & bsButton
observe({
    x <- input$box_globalCutoff
    toggle(id="globalCutoff",  condition=x)
    toggle(id="button_globalCutoff", condition=x)
})

# renderUI if checkboxInput == 1
output$globalCutoffUI <- renderUI({
    globalCutoffUI
})

# set global cutoff upon bsButton click
observeEvent(input$button_globalCutoff, {
    sep_cutoffs(vals$dbFrame2) <- as.numeric(input$input_globalCutoff)
})

# sliderInput: "Mahalanobis distance threshold"
mahalCutoff <- reactive(input$mahalCutoff)

# apply cutoffs if deconvolution parameters change
dbFrame <- reactive({
    x <- vals$dbFrame2
    if (is.null(x)) return()
    applyCutoffs(x, mahalCutoff())
})

# ------------------------------------------------------------------------------
# yieldPlot, eventPlot & mahalPlot
# ------------------------------------------------------------------------------

# inputSelect choices
yieldPlotChoices <- reactive({
    x <- dbFrame()
    if (is.null(x)) {
        NULL
    } else {
        ids <- rownames(debaKey())
        setNames(c(0, ids), c("All", ids))
    }
})
eventPlotChoices <- reactive({
    x <- dbFrame()
    if (is.null(x)) {
        NULL
    } else {
        sort(unique(bc_ids(dbFrame())))
    }
})
mahalPlotChoices <- reactive({
    x <- dbFrame()
    if (is.null(x)) {
        NULL
    } else {
        choices <- eventPlotChoices()
        choices[choices != 0]
    }
})

# render IUs
output$yieldPlotPanel <- renderUI({
    x <- yieldPlotChoices()
    if (!is.null(x)) 
        yieldPlotPanel(x)
})
output$eventPlotPanel <- renderUI({
    x <- eventPlotChoices()
    if (!is.null(x)) 
        eventPlotPanel(x)
})
output$mahalPlotPanel <- renderUI({
    x <- mahalPlotChoices()
    if (!is.null(x)) 
        mahalPlotPanel(x)
})

# get selected ID
yieldPlotSelected <- reactive(input$yieldPlotSelected)
eventPlotSelected <- reactive(input$eventPlotSelected)
mahalPlotSelected <- reactive(input$mahalPlotSelected)

# get n_events & cofactor for eventPlot & mahalPlot
eventPlotNEvents  <- reactive(as.numeric(input$n_events))
mahalPlotCofactor <- reactive(input$mahalPlotCofactor)

# render plots
output$yieldPlot <- renderPlot({
    if (is.null(dbFrame())) 
        return()
    plotYields(
        dbFrame(), 
        yieldPlotSelected())
})
output$eventPlot <- renderPlot({
    if (is.null(dbFrame()))
        return()
    plotEvents(
        dbFrame(), 
        eventPlotSelected(), 
        eventPlotNEvents())
})
output$mahalPlot <- renderPlot({
    if (is.null(dbFrame()))
        return()
    plotMahal(
        dbFrame(), 
        mahalPlotSelected(), 
        mahalPlotCofactor())
})

# renderDataTable: IDs | Counts | Cutoffs | Yields
output$table_summary <- DT::renderDataTable({ 
    if (is.null(dbFrame())) return()
    summary_tbl(dbFrame()) 
})

# ------------------------------------------------------------------------------
# next / previous buttons
# ------------------------------------------------------------------------------
observeEvent(input$prev_yieldPlot, { 
    x <- yieldPlotChoices()
    y <- which(x == input$select_yieldPlot)
    if (y == 1) return()
    updateSelectInput(session, "select_yieldPlot", selected=x[y-1]) 
})
observeEvent(input$next_yieldPlot, { 
    x <- yieldPlotChoices()
    y <- which(x == input$select_yieldPlot)
    if (y == length(x)) return()
    updateSelectInput(session, "select_yieldPlot", selected=x[y+1]) 
})

observeEvent(input$prev_eventPlot, { 
    x <- eventPlotChoices()
    y <- which(x == input$select_eventPlot)
    if (y == 1) return()
    updateSelectInput(session, "select_eventPlot", selected=x[y-1]) 
})
observeEvent(input$next_eventPlot, { 
    x <- eventPlotChoices()
    y <- which(x == input$select_eventPlot)
    if (y == length(x)) return()
    updateSelectInput(session, "select_eventPlot", selected=x[y+1])
})

observeEvent(input$prev_mahalPlot, { 
    x <- mahalPlotChoices()
    y <- which(x == input$select_mahalPlot)
    if (y == 1) return()
    updateSelectInput(session, "select_mahalPlot", selected=x[y-1])
})
observeEvent(input$next_mahalPlot, { 
    x <- mahalPlotChoices()
    y <- which(x == input$select_mahalPlot)
    if (y == length(x)) return()
    updateSelectInput(session, "select_mahalPlot", selected=x[y+1])
})

# ------------------------------------------------------------------------------
# synchronize eventPlot yieldPlot and mahalPlot
# ------------------------------------------------------------------------------
observe({
    x <- input$select_yieldPlot
    if (is.null(x) || x == 0) return() 
    updateSelectInput(session, "select_eventPlot", selected=x)
    updateSelectInput(session, "select_mahalPlot", selected=x)
})
observe({
    x <- input$select_eventPlot
    if (is.null(x) || x == 0) return() 
    updateSelectInput(session, "select_yieldPlot", selected=x)
    updateSelectInput(session, "select_mahalPlot", selected=x)
})
observe({
    x <- input$select_mahalPlot
    if (is.null(x)) return()
    updateSelectInput(session, "select_yieldPlot", selected=x)
    updateSelectInput(session, "select_eventPlot", selected=x)
})
