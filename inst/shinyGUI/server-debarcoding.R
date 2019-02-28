# ==============================================================================
# DEBARCODING
# ==============================================================================

# read input FCS
ffDeba <- reactive({
    if (vals$keepDataComp) {
        fs <- fsComped()
        print(fs)
        
        CATALYST::concatFCS(fs)
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
    match <- grep("(.csv)$", input$debaSchemeCsv$name, ignore.case=TRUE)
    isCSV <- length(match) == 1
    if (!isCSV) {
        showNotification(
            h4(strong("Input debarcoding scheme should be a CSV file.")),
            duration=NULL, type="error")
        return()
    }
    key <- read.csv(input$debaSchemeCsv$datapath, 
        check.names=FALSE, row.names=1)
    vals$debaKeyIsValid <- check_key(key, ffDeba())
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
    disable(id="debarcodeDeba")
    # assignPrelim()
    showNotification(h4(strong("Assigning preliminary IDs...")), 
        duration=NULL, closeButton=FALSE, id="msg", type="message")
    vals$dbFrame1Deba <- CATALYST::assignPrelim(ffDeba(), debaKey())
    removeNotification(id="msg")
    # estCutoffs()
    showNotification("Estimating separation cutoffs...",
        id="estimating_sep_cutoffs", duration=NULL, closeButton=FALSE)
    vals$cutoff_ests_deba <- sep_cutoffs(CATALYST::estCutoffs(x=vals$dbFrame1Deba))
    removeNotification(id="estimating_sep_cutoffs")
    # extend sidebar
    output$debarcodingSidebar2 <- renderUI(debarcodingSidebar2)
    enable(id="debarcodeDeba")
})

# ------------------------------------------------------------------------------
# cutoff adjustment
# ------------------------------------------------------------------------------
# render UI for cutoff adjustment or input for global cutoff
output$deba_cutoffs_UIDeba <- renderUI({
    switch(input$deba_cutoffsDeba,
        est_cutoffs=NULL,
        adj_cutoffs=
            adjustCutoffUI(
                dbFrame=vals$dbFrame2Deba, 
                choices=adjustCutoffChoicesDeba(),
                module="Deba"),
        global_cutoff=
            globalCutoffUI(module="Deba"))
})

# use cutoff estimates
observeEvent(input$deba_cutoffsDeba == "est_cutoffs", ignoreInit=TRUE, {
    sep_cutoffs(vals$dbFrame1Deba) <- vals$cutoff_ests_deba
    vals$dbFrame2Deba <- vals$dbFrame1Deba
})

# get selectInput choices for cutoff adjustment
adjustCutoffChoicesDeba <- reactive({
    req(dbFrameDeba())
    rownames(bc_key(dbFrameDeba()))
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

# ------------------------------------------------------------------------------
# download handlers
# ------------------------------------------------------------------------------

smplNmsDeba <- reactive({
    req(dbFrameDeba())
    if (input$box_IDsAsNms) {
        c("Unassigned", rownames(bc_key(dbFrameDeba())))
    } else if (input$box_upldNms) {
        req(input$input_upldNms)
        read.csv(input$input_upldNms$datapath, header=FALSE)
    }
})

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
        tmpdir <- tempdir()
        setwd(tmpdir)
        # ----------------------------------------------------------------------
        # debarcoding summary table: 
        # IDs | Counts | Cutoffs | Yields
        ids <- rownames(bc_key(dbFrameDeba()))
        nBcs <- nrow(bc_key(dbFrameDeba()))
        cutoffs <- sep_cutoffs(dbFrameDeba())
        yields <- yields(dbFrameDeba())[cbind(seq_len(nBcs), 
            findInterval(cutoffs, seq(0, 1, .01)))]
        yields <- round(100*yields, 4)
        tbl <- matrix(cbind(ids, sapply(ids, function(id) 
            sum(bc_ids(dbFrameDeba()) == id)), cutoffs, yields), ncol=4,
            dimnames=list(NULL, c("ID", "Count","Cutoff", "Yield [%]")))
        tblNm <- paste0(format(Sys.Date(), "%y%m%d"), "-debarcoding.csv") 
        write.csv(tbl, file.path(tmpdir, tblNm), row.names=FALSE)
        # ----------------------------------------------------------------------
        # population-wise FCS files
        # ----------------------------------------------------------------------
        # get output file names
        if (input$box_IDsAsNms) {
            smplNms <- smplNmsDeba()
        } else if (input$box_upldNms) {
            # check that a name has been supplied for all sample
            # & that sample IDs match with the input barcoding scheme
            if (nrow(smplNmsDeba()) < nBcs) {
                showNotification(paste("Only", nrow(smplNmsDeba()), 
                    "sample names provided but", nBcs, "needed."),
                    type="error", closeButton=FALSE)
                return()
            } else if (sum(smplNmsDeba()[, 1] %in% ids) != nBcs) {
                showNotification(
                    "Couldn't find a file name for all samples.\n
                    Please make sure all sample IDs occur\n
                    in the provided naming scheme.",
                    type="error", closeButton=FALSE)
                return()
            }
            smplNms <- c("Unassigned", paste0(smplNmsDeba()[, 2], "_", ids))
        }
        unique_ids <- unique(bc_ids(dbFrameDeba()))
        nFiles <- length(unique_ids)
        inds <- match(bc_ids(dbFrameDeba()), unique_ids)
        out_nms <- paste0(smplNms, ".fcs")
        out_nms <- out_nms[match(unique_ids, c(0, ids))]
        # match assignments with IDs
        # write population-wise FCS file
        withProgress(message="Debarcoding samples...", value=0, {
            for (i in seq_along(unique_ids)) {
                suppressWarnings(flowCore::write.FCS(
                    x=ffDeba()[inds == i, ], 
                    filename=file.path(tmpdir, out_nms[i])))
                incProgress(1/nFiles, detail=paste0(i, "/", nFiles))
            }
        })
        showNotification(h4(strong("Writing FCS files...")),
            id="msg", duration=NULL, closeButton=NULL, type="default")
        fileNms <- c(tblNm, smplNms)
        zip(zipfile=file, files=out_nms) 
        removeNotification(id="msg")
        }, 
    contentType="application/zip")                        

output$dwnld_debaPlots <- downloadHandler(
    filename=function() { "yield_event_plots.zip" },
    content =function(file) { 
        tmpdir <- tempdir()
        setwd(tmpdir)
        CATALYST::plotYields(
            x=dbFrameDeba(), 
            which=yieldPlotChoices(), 
            out_path=tmpdir)
        CATALYST::plotEvents(
            x=dbFrameDeba(), 
            out_path=tmpdir, 
            n_events=250)
        zip(zipfile=file,
            files=paste0(c("yield_plot", "event_plot"), ".pdf")) },
    contentType="application/zip")
