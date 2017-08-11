# ==============================================================================
# DEBARCODING
# ==============================================================================

# read input FCS
ffDeba <- reactive({
    if (vals$keepDataComp) {
        ff <- fsComped()
        if (input$box_setToZero)
            exprs(ff)[exprs(ff) < 0] <- 0
        ff
    } else {
        req(input$fcsDeba)
        # check validity of input FCS files
        valid <- check_FCS_fileInput(input$fcsDeba)
        if (!valid) return()
        read.FCS(
            filename=input$fcsDeba$datapath,
            transformation=FALSE,
            truncate_max_range=FALSE)
    }
})


# expand sidebar
output$debarcodingSidebar1 <- renderUI({
    req(ffDeba())
    debarcodingSidebar1
})

# fileInput "upload barcoding scheme (CSV)"
output$selectBcChs <- renderUI({
    req(input$boxSelectBcChs == 1)
    chs <- flowCore::colnames(ffDeba())
    ms <- as.numeric(gsub("[[:alpha:][:punct:]]", "", chs))
    chs <- chs[!is.na(as.numeric(ms))]
    selectInput(
        inputId="input_bcChs", label=NULL, choices=chs, 
        multiple=TRUE, selectize=FALSE, size=12)
})

# check validity of debarcoding scheme CSV
observe({
    req(input$debaSchemeCsv)
    key <- read.csv(input$debaSchemeCsv$datapath, 
        check.names=FALSE, row.names=1)
    vals$debaKeyIsValid <- checkKey(key, ffDeba())
})

# get debarcoding scheme
debaKey <- reactive({
    req(isTRUE(vals$debaKeyIsValid))
    key <- read.csv(input$debaSchemeCsv$datapath, 
        check.names=FALSE, row.names=1)
})

# ------------------------------------------------------------------------------
# bsButton: "Debarcode"
# ------------------------------------------------------------------------------
observe({
    toggleState(id="debarcodeDeba", condition=req(isTRUE(vals$debaKeyIsValid)))
})

observeEvent(input$debarcodeDeba, {
    showNotification(h5("Assigning preliminary IDs..."), 
        duration=NULL, closeButton=FALSE, id="msg", type="message")
    vals$dbFrame1Deba <- CATALYST::assignPrelim(ffDeba(), debaKey())
    removeNotification(id="msg")
    output$debarcodingSidebar2 <- renderUI(debarcodingSidebar2)
})

# ------------------------------------------------------------------------------
# toggle checkboxes
observe({
    req(input$checkbox_estCutoffsDeba == 1)
    updateCheckboxInput(session, "checkbox_adjustCutoffDeba", value=FALSE)
    updateCheckboxInput(session, "checkbox_globalCutoffDeba", value=FALSE)
})
observe({
    req(input$checkbox_adjustCutoffDeba == 1)
    updateCheckboxInput(session, "checkbox_estCutoffsDeba",   value=FALSE)
    updateCheckboxInput(session, "checkbox_globalCutoffDeba", value=FALSE)
})
observe({
    req(input$checkbox_globalCutoffDeba == 1)
    updateCheckboxInput(session, "checkbox_estCutoffsDeba",   value=FALSE)
    updateCheckboxInput(session, "checkbox_adjustCutoffDeba", value=FALSE)
})

# estCutoffs()
observeEvent(input$checkbox_estCutoffsDeba, {
    vals$dbFrame2Deba <- CATALYST::estCutoffs(x=vals$dbFrame1Deba)
})

# ------------------------------------------------------------------------------
# cutoff adjustment
# ------------------------------------------------------------------------------

# toggle adjustCutoffsUI
observe({
    toggle(
        id="adjustCutoffUIDeba", 
        condition=input$checkbox_adjustCutoffDeba)
})

# get selectInput choices
adjustCutoffChoicesDeba <- reactive({
    req(dbFrameDeba())
    rownames(bc_key(dbFrameDeba()))
})

# render adjustCutoffsUI
output$adjustCutoffUIDeba <- renderUI({
    adjustCutoffUI(
        dbFrame=vals$dbFrame2Deba, 
        choices=adjustCutoffChoicesDeba(),
        module="Deba")
})

# synchronize selectInput & numericInput w/ yield plot
observe({
    x <- selectedYieldPlotDeba()
    req(!is.null(x), x != 0)
    updateSelectInput(session, "select_adjustCutoffDeba", selected=x)
    updateNumericInput(session, "input_adjustCutoffDeba", 
        value=paste(sep_cutoffs(vals$dbFrame2Deba)[x]))
})
observe({
    req(input$select_adjustCutoffDeba)
    updateSelectInput(session, "select_yieldPlotDeba", 
        selected=input$select_adjustCutoffDeba)
})

# adjust cutoff upon bsButton click
observeEvent(input$button_adjustCutoffDeba, {
    x <- match(input$select_adjustCutoffDeba, adjustCutoffChoicesDeba())
    sep_cutoffs(vals$dbFrame2Deba)[x] <- 
        as.numeric(input$input_adjustCutoffDeba)
})

# ------------------------------------------------------------------------------
# global separation cutoff
# ------------------------------------------------------------------------------

# toggle globalCutoffUI
observe({
    toggle(
        id="globalCutoffUIDeba", 
        condition=input$checkbox_globalCutoffDeba)
})

# render globalCutoffUI
output$globalCutoffUIDeba <- renderUI(globalCutoffUI(module="Deba"))

# set global cutoff upon bsButton click
observeEvent(input$button_globalCutoffDeba, {
    sep_cutoffs(vals$dbFrame2Deba) <- as.numeric(input$input_globalCutoffDeba)
})

# sliderInput: "Mahalanobis distance threshold"
observeEvent(input$button_mhlCutoffDeba, {
    vals$mhlCutoffDeba <- input$mhlCutoffDeba
})

# apply cutoffs if deconvolution parameters change
dbFrameDeba <- reactive({
    req(vals$dbFrame2Deba)
    CATALYST::applyCutoffs(x=vals$dbFrame2Deba, mhl_cutoff=vals$mhlCutoffDeba)
})

# ------------------------------------------------------------------------------
# yieldPlot, eventPlot & mahalPlot
# ------------------------------------------------------------------------------

# inputSelect choices
yieldPlotChoices <- reactive({
    req(vals$dbFrame1Deba)
    ids <- rownames(bc_key(vals$dbFrame1Deba))
    setNames(c(0, ids), c("All", ids))
})
eventPlotChoices <- reactive({
    req(vals$dbFrame2Deba)
    sort(unique(bc_ids(vals$dbFrame2Deba)))
})
mahalPlotChoices <- reactive({
    req(vals$dbFrame2Deba)
    choices <- eventPlotChoices()
    choices[choices != 0]
})

# render UIs
output$yieldPlotPanelDeba <- renderUI({
    req(yieldPlotChoices())
    yieldPlotModule(yieldPlotChoices(), "Deba")
})
output$eventPlotPanel <- renderUI({
    x <- eventPlotChoices()
    y <- selectedYieldPlotDeba()
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
output$yieldPlotDeba <- renderPlotly({
    req(dbFrameDeba())
    plotYields(
        x=dbFrameDeba(), 
        which=selectedYieldPlotDeba())[[1]]
})
output$eventPlot <- renderPlot({
    req(dbFrameDeba())
    plotEvents(
        x=dbFrameDeba(), 
        which=input$select_eventPlot, 
        n_events=eventPlotNEvents())
})
output$mahalPlot <- renderPlot({
    req(dbFrameDeba())
    plotMahal(
        x=dbFrameDeba(), 
        which=input$select_mahalPlot, 
        cofactor=mahalPlotCofactor())
})

# renderDataTable: IDs | Counts | Cutoffs | Yields
output$summaryTblDeba <- DT::renderDataTable({ 
    req(dbFrameDeba())
    summary_tbl(dbFrameDeba()) 
})

# keep track of currently selected sample
selectedYieldPlotDeba <- reactive({input$select_yieldPlotDeba})

# ------------------------------------------------------------------------------
# next / previous buttons
# ------------------------------------------------------------------------------
observe({
    choices <- yieldPlotChoices()
    n <- length(choices)
    selected <- match(selectedYieldPlotDeba(), choices)
    toggleState(id="prev_yieldPlotDeba", condition=selected != 1)
    toggleState(id="next_yieldPlotDeba", condition=selected != n)
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

observeEvent(input$prev_yieldPlotDeba, { 
    choices <- yieldPlotChoices()
    selected <- match(selectedYieldPlotDeba(), choices)
    updateSelectInput(session, 
        inputId="select_yieldPlotDeba", 
        selected=choices[selected-1]) 
})
observeEvent(input$next_yieldPlotDeba, { 
    choices <- yieldPlotChoices()
    selected <- match(selectedYieldPlotDeba(), choices)
    updateSelectInput(session, 
        inputId="select_yieldPlotDeba", 
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
    x <- selectedYieldPlotDeba()
    if (is.null(x) || x == 0) return() 
    updateSelectInput(session, "select_eventPlot", selected=x)
    updateSelectInput(session, "select_mahalPlot", selected=x)
})
observe({
    x <- input$select_eventPlot
    if (is.null(x) || x == 0) return() 
    updateSelectInput(session, "select_yieldPlotDeba", selected=x)
    updateSelectInput(session, "select_mahalPlot", selected=x)
})
observe({
    x <- input$select_mahalPlot
    if (is.null(x)) return()
    updateSelectInput(session, "select_yieldPlotDeba", selected=x)
    updateSelectInput(session, "select_eventPlot", selected=x)
})

# bsButton "Go to compensation"
observeEvent(input$goToComp, {
    shinyjs::hide(id="fcsComp")
    updateTabItems(session, inputId="tabs", selected="compensation")
    vals$keepDataDeba <- TRUE
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

smplNmsDeba <- reactive({
    req(dbFrameDeba())
    if (input$box_IDsAsNms) {
        c("Unassigned", rownames(bc_key(dbFrameDeba())))
    } else if (input$box_upldNms) {
        req(input$input_upldNms)
        read.csv(input$input_upldNms$datapath, header=FALSE)
    }
})

# toggle fileInput: "Upload naming sheet (CSV)"
observe(toggle(id="input_upldNms", condition=input$box_upldNms))

# toggle downloadButton
observe({
    req(dbFrameDeba())
    test <- input$box_IDsAsNms || (
        input$box_upldNms == 1 && !is.null(input$input_upldNms))
    toggleState(id="dwnld_debaFCS", condition=test)
})

output$dwnld_debaFcs <- downloadHandler(
    filename=function() { 
        paste0(format(Sys.Date(), "%y%m%d"), "-debarcoding.zip")
    },
    content =function(file) { 
        nms <- smplNmsDeba()
        dbFrame <- dbFrameDeba()
        ids <- rownames(bc_key(dbFrame))
        inds <- c(0, ids) %in% sort(unique(bc_ids(dbFrame)))
        key <- debaKey()
        tmpdir <- tempdir()
        setwd(tmpdir)
        
        if (input$box_IDsAsNms) {
            CATALYST::outFCS(
                x=dbFrame, 
                y=ffDeba(), 
                out_path=tmpdir) 
            fileNms <- nms[inds]
        } else if (input$box_upldNms) {
            if (nrow(nms) < nrow(key)) {
                showNotification(paste("Only", nrow(nms), 
                    "sample names provided but", nrow(key), "needed."),
                    type="error", closeButton=FALSE)
                return()
            } else if (sum(nms[, 1] %in% ids) != length(ids)) {
                showNotification(
                    "Couldn't find a file name for all samples.\n
                    Please make sure all sample IDs occur\n
                    in the provided naming scheme.",
                    type="error", closeButton=FALSE)
                return()
            }
            CATALYST::outFCS(
                x=dbFrame, 
                y=ffDeba(),
                out_path=tmpdir, 
                out_nms=paste0(nms[, 2], "_", ids))
            fileNms <- c("Unassigned", paste0(nms[, 2], "_", ids))[inds]
        }
        write.csv()
        fileNms <- c("summary_table.csv", paste0(fileNms, ".fcs"))
        zip(zipfile=file, files=fileNms) 
        }, 
    contentType="application/zip")                        

output$dwnld_debaPlots <- downloadHandler(
    filename=function() { "yield_event_plots.zip" },
    content =function(file) { 
        dbFrame <- dbFrameDeba
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