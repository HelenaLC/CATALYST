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
bcChoices <- reactive({

})
output$selectBcChs <- renderUI({
    req(input$boxSelectBcChs == 1)
    chs <- flowCore::colnames(ffDeba())
    ms <- as.numeric(gsub("[[:alpha:][:punct:]]", "", chs))
    chs <- chs[!is.na(as.numeric(ms))]
    selectInput(
        inputId="input_bcChs", label=NULL, choices=chs, 
        multiple=TRUE, selectize=FALSE, size=12)
})

# get debarcoding scheme
debaKey <- reactive({
    if (!is.null(input$input_bcChs)) {
        as.numeric(gsub("[[:alpha:][:punct:]]", "", input$input_bcChs))
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

# ------------------------------------------------------------------------------
# cutoff adjustment
# ------------------------------------------------------------------------------

# toggle UI
observe({
    x <- input$box_adjustCutoff
    toggle(id="adjustCutoffUI", condition=x)
})

# selectInput choices
adjustCutoffChoices <- reactive({
    req(dbFrame())
    rownames(bc_key(dbFrame()))
})

# renderUI if checkboxInput == 1
output$adjustCutoffUI <- renderUI({
    adjustCutoffUI(
        dbFrame=vals$dbFrame2, 
        choices=adjustCutoffChoices())
})

# synchronize selectInput & numericInput w/ yield plot
observe({
    x <- input$select_yieldPlot
    req(!is.null(x), x != 0)
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

# toggle UI
observe({
    x <- input$box_globalCutoff
    toggle(id="globalCutoffUI", condition=x)
})

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
    req(vals$dbFrame1)
    ids <- rownames(bc_key(vals$dbFrame1))
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
    req(x)
    yieldPlotPanel(x)
})
output$eventPlotPanel <- renderUI({
    x <- eventPlotChoices()
    y <- input$select_yieldPlot
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
        x=dbFrame(), 
        which=input$select_yieldPlot)
})
output$eventPlot <- renderPlot({
    req(dbFrame())
    plotEvents(
        x=dbFrame(), 
        which=input$select_eventPlot, 
        n_events=eventPlotNEvents())
})
output$mahalPlot <- renderPlot({
    req(dbFrame())
    plotMahal(
        x=dbFrame(), 
        which=input$select_mahalPlot, 
        cofactor=mahalPlotCofactor())
})

# renderDataTable: IDs | Counts | Cutoffs | Yields
output$table_summary <- DT::renderDataTable({ 
    req(dbFrame())
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

# ------------------------------------------------------------------------------
# download handlers
# ------------------------------------------------------------------------------

# toggle checkboxes
observe({
    req(input$box_IDsAsNms == 1)
    updateCheckboxInput(session, "box_upldNms", value=FALSE)
})
observe({
    req(input$box_upldNms == 1)
    updateCheckboxInput(session, "box_IDsAsNms", value=FALSE)
})

# toggle fileInput: "Upload naming sheet (CSV)"
observe(toggle(id="input_upldNms", condition=input$box_upldNms))

output$dwnld_debaFcs <- downloadHandler(
    filename=function() { "fcs.zip" },
    content =function(file) { 
        dbFrame <- dbFrame()
        key <- debaKey()
        ids <- rownames(bc_key(dbFrame))
        tmpdir <- tempdir()
        setwd(tmpdir)
        inds <- c(0, ids) %in% sort(unique(bc_ids(dbFrame)))
        if (input$box_IDsAsNms) {
            CATALYST::outFCS(
                x=dbFrame, 
                y=ffDeba(), 
                out_path=tmpdir) 
            fileNms <- c("Unassigned", ids)[inds]
        } else if (input$box_upldNms) {
            if (is.null(input$input_upldNms)) {
                showNotification("Please upload a naming.",
                    type="error", closeButton=FALSE)
                return()
            } 
            smplNms <- read.csv(input$input_upldNms$datapath, header=FALSE)
            if (nrow(smplNms) < nrow(key)) {
                showNotification(paste("Only", nrow(smplNms), 
                    "sample names provided but", nrow(key), "needed."),
                    type="error", closeButton=FALSE)
                return()
            } else if (sum(smplNms[, 1] %in% ids) != nrow(key)) {
                showNotification(
                    "Couldn't find a file name for all samples.\n
                    Please make sure all sample IDs occur\n
                    in the provided naming scheme.",
                    type="error", closeButton=FALSE)
                return()
            }
            CATALYST::outFCS(
                x=dbFrame(), 
                y=ffDeba(),
                out_path=tmpdir, 
                out_nms=paste0(smplNms[, 2], "_", ids))
            fileNms <- c("Unassigned", paste0(smplNms[, 2], "_", ids))[inds]
        }
        zip(zipfile=file, 
            files=paste0(fileNms, ".fcs")) }, 
    contentType="application/zip")                        

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