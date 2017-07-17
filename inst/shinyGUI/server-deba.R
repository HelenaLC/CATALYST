# ==============================================================================
# DEBARCODING
# ==============================================================================

# read input FCS
ffDeba <- reactive(flowCore::read.FCS(input$fcsDeba$datapath))

# expand sidebar
output$debarcodingSidebar1 <- renderUI({
    req(input$fcsDeba)
    debarcodingSidebar1
})

# ------------------------------------------------------------------------------
# checkboxInputs 
# "Upload barcoding scheme (CSV)" &
# "Select single-positive channels"
# ------------------------------------------------------------------------------

# toggle checkboxes
observe({
    req(input$boxUploadCsv == 1)
    updateCheckboxInput(session, "boxSelectBcChs", value=FALSE)
})
observe({
    req(input$boxSelectBcChs == 1)
    updateCheckboxInput(session, "boxUploadCsv", value=FALSE)
})

# boxUploadCsv -> render fileInput
output$debaSchemeCsv <- renderUI({
    req(input$boxUploadCsv == 1)
    fileInput(inputId="debaSchemeCsv", label=NULL, accept=".csv")
})

# boxSelectBcChs -> render selectInput
output$selectBcChs <- renderUI({
    req(input$boxSelectBcChs == 1)
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
    }
})

# ------------------------------------------------------------------------------
# bsButton: "Debarcode"
# ------------------------------------------------------------------------------
observe({
    toggleState(id="buttonDebarcode", condition=!is.null(debaKey()))
})

observeEvent(input$buttonDebarcode, {
    showNotification(h5("Assigning preliminary IDs..."), 
        duration=NULL, closeButton=FALSE, id="msg", type="message")
    vals$dbFrame1 <- assignPrelim(ffDeba(), debaKey())
    removeNotification(id="msg")
    output$debarcodingSidebar2 <- renderUI(debarcodingSidebar2)
})

# ------------------------------------------------------------------------------
# toggle checkboxes
observe({
    req(input$box_estCutoffs == 1)
    updateCheckboxInput(session, "box_adjustCutoff", value=FALSE)
    updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
})
observe({
    req(input$box_adjustCutoff == 1)
    updateCheckboxInput(session, "box_estCutoffs",   value=FALSE)
    updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
})
observe({
    req(input$box_globalCutoff == 1)
    updateCheckboxInput(session, "box_estCutoffs",   value=FALSE)
    updateCheckboxInput(session, "box_adjustCutoff", value=FALSE)
})

# ------------------------------------------------------------------------------
# cutoff estimation
# ------------------------------------------------------------------------------
observeEvent(input$box_estCutoffs, {
    vals$dbFrame2 <- estCutoffs(vals$dbFrame1)
})

# toggle UIs
observe({
    x <- input$box_adjustCutoff
    toggle(id="adjustCutoffUI", condition=x)
})
observe({
    x <- input$box_globalCutoff
    toggle(id="globalCutoffUI", condition=x)
})

# ------------------------------------------------------------------------------
# cutoff adjustment
# ------------------------------------------------------------------------------

# selectInput choices
adjustCutoffChoices <- reactive(rownames(debaKey()))

# renderUI if checkboxInput == 1
output$adjustCutoffUI <- renderUI({
    adjustCutoffUI(
        dbFrame=vals$dbFrame2, 
        choices=adjustCutoffChoices(), 
        selected=1)
})

# synchronize selectInput & numericInput w/ yield plot
observe({
    x <- input$select_yieldPlot
    req(x != 0)
    updateSelectInput(session, "select_adjustCutoff", selected=x)
    updateNumericInput(session, "input_adjustCutoff", 
        value=paste(sep_cutoffs(vals$dbFrame2)[x]))
})
observe({
    req(input$select_adjustCutoff)
    updateSelectInput(session, "select_yieldPlot", 
        selected=input$select_adjustCutoff)
})

# adjust cutoff upon bsButton click
observeEvent(input$button_adjustCutoff, {
    x <- match(input$select_adjustCutoff, adjustCutoffChoices())
    sep_cutoffs(vals$dbFrame2)[x] <- as.numeric(input$input_adjustCutoff)
})

# ------------------------------------------------------------------------------
# global separation cutoff
# ------------------------------------------------------------------------------

# renderUI if checkboxInput == 1
output$globalCutoffUI <- renderUI(globalCutoffUI)

# set global cutoff upon bsButton click
observeEvent(input$button_globalCutoff, {
    sep_cutoffs(vals$dbFrame2) <- as.numeric(input$input_globalCutoff)
})

# sliderInput: "Mahalanobis distance threshold"
observeEvent(input$button_mhlCutoff, {
    vals$mhlCutoffDeba <- input$input$mahalCutoff
})

# apply cutoffs if deconvolution parameters change
dbFrame <- reactive({
    req(vals$dbFrame2)
    applyCutoffs(vals$dbFrame2, vals$mhlCutoffDeba)
})

# ------------------------------------------------------------------------------
# yieldPlot, eventPlot & mahalPlot
# ------------------------------------------------------------------------------

# inputSelect choices
yieldPlotChoices <- reactive({
    req(debaKey())
    ids <- rownames(debaKey())
    setNames(c(0, ids), c("All", ids))
})
eventPlotChoices <- reactive({
    req(vals$dbFrame2)
    sort(unique(bc_ids(vals$dbFrame2)))
})
mahalPlotChoices <- reactive({
    req(vals$dbFrame2)
    choices <- eventPlotChoices()
    choices[choices != 0]
})

# render UIs
output$yieldPlotPanel <- renderUI({
    x <- yieldPlotChoices()
    y <- input$select_yieldPlot
    if (!is.null(x)) 
        yieldPlotPanel(x, y)
})
output$eventPlotPanel <- renderUI({
    x <- eventPlotChoices()
    y <- input$select_eventPlot
    if (!is.null(x)) 
        eventPlotPanel(x, y)
})
output$mahalPlotPanel <- renderUI({
    x <- mahalPlotChoices()
    y <- input$select_mahalPlot
    if (!is.null(x)) 
        mahalPlotPanel(x, y)
})

# get n_events & cofactor for eventPlot & mahalPlot
eventPlotNEvents  <- reactive(as.numeric(input$n_events))
mahalPlotCofactor <- reactive(input$mahalPlotCofactor)

# render plots
output$yieldPlot <- renderPlot({
    req(dbFrame())
    plotYields(
        dbFrame(), 
        input$select_yieldPlot)
})
output$eventPlot <- renderPlot({
    req(dbFrame())
    plotEvents(
        dbFrame(), 
        input$select_eventPlot, 
        eventPlotNEvents())
})
output$mahalPlot <- renderPlot({
    req(dbFrame())
    plotMahal(
        dbFrame(), 
        input$select_mahalPlot, 
        mahalPlotCofactor())
})

# renderDataTable: IDs | Counts | Cutoffs | Yields
output$table_summary <- DT::renderDataTable({ 
    if (!is.null(dbFrame()))
        summary_tbl(dbFrame()) 
})

# ------------------------------------------------------------------------------
# next / previous buttons
# ------------------------------------------------------------------------------
observe({
    choices <- yieldPlotChoices()
    n <- length(choices)
    selected <- match(input$select_yieldPlot, choices)
    toggleState(id="prev_yieldPlot", condition=selected != 1)
    toggleState(id="next_yieldPlot", condition=selected != n)
})
observe({
    choices <- eventPlotChoices()
    n <- length(choices)
    selected <- match(input$select_eventPlot, choices)
    toggleState(id="prev_eventPlot", condition=selected != 1)
    toggleState(id="next_eventPlot", condition=selected != n)
})
observe({
    choices <- mahalPlotChoices()
    n <- length(choices)
    selected <- match(input$select_mahalPlot, choices)
    toggleState(id="prev_mahalPlot", condition=selected != 1)
    toggleState(id="next_mahalPlot", condition=selected != n)
})

observeEvent(input$prev_yieldPlot, { 
    choices <- yieldPlotChoices()
    selected <- match(input$select_yieldPlot, choices)
    updateSelectInput(session, 
        inputId="select_yieldPlot", 
        selected=choices[selected-1]) 
})
observeEvent(input$next_yieldPlot, { 
    choices <- yieldPlotChoices()
    selected <- match(input$select_yieldPlot, choices)
    updateSelectInput(session, 
        inputId="select_yieldPlot", 
        selected=choices[selected+1]) 
})

observeEvent(input$prev_eventPlot, { 
    choices <- eventPlotChoices()
    selected <- match(input$select_eventPlot, choices)
    updateSelectInput(session, 
        inputId="select_eventPlot", 
        selected=choices[selected-1]) 
})
observeEvent(input$next_eventPlot, { 
    choices <- eventPlotChoices()
    selected <- match(input$select_eventPlot, choices)
    updateSelectInput(session, 
        inputId="select_eventPlot", 
        selected=choices[selected+1])
})

observeEvent(input$prev_mahalPlot, { 
    choices <- mahalPlotChoices()
    selected <- match(input$select_mahalPlot, choices)
    updateSelectInput(session, 
        inputId="select_mahalPlot", 
        selected=choices[selected-1])
})
observeEvent(input$next_mahalPlot, { 
    choices <- mahalPlotChoices()
    selected <- match(input$select_mahalPlot, choices)
    updateSelectInput(session, 
        inputId="select_mahalPlot", 
        selected=choices[selected+1])
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
