# ==============================================================================
# COMPENSATION
# ==============================================================================

# get flowSet to compensated: 
# initially, use normalized samples or read input FCS file(s)
# once duplicate masses or metal mismatches where checked,
# the newly constructed data will be used
fsComp <- reactive({
    data <- NULL
    if (!is.null(fs <- vals$fsComp_metsChecked)) {
        data <- fs
    } else if (!is.null(fs <- vals$fsComp_msChecked)) {
        data <- fs
    } else if (vals$keepDataNorm) {
        n <- length(ffsNormed())
        if (n == 1) {
            ff <- ffsNormed()[[1]]
        } else {
            ff <- CATALYST::concatFCS(x=ffsNormed())
        }
        if (input$box_removeBeads) {
            removed <- sapply(seq_len(n), function(i) 
                mhlDists()[[i]] < vals$mhlCutoffs[i])
            data <- list(ff[!removed, ])
        } else {
            data <- list(ff)
        }
    } else {
        req(input$fcsComp)
        n <- nrow(input$fcsComp)
        # check validity of input FCS files
        valid <- check_FCS_fileInput(input$fcsComp, n)
        if (!valid) return()
        data <- lapply(seq_len(n), function(i) 
            flowCore::read.FCS(
                filename=input$fcsComp[[i, "datapath"]],
                transformation=FALSE,
                truncate_max_range=FALSE))
    }
    req(data)
    as(data, "flowSet")
})

# render fileInput for spillover matrix CSV
output$uploadSpillMat <- renderUI({
    req(input$checkbox_upldSm == 1)
    fileInput(
        inputId="inputSpillMat", 
        label=NULL, 
        accept=".csv")
})

# render fileInput for single-stained controls
output$uploadSs <- renderUI({
    req(input$checkbox_estSm == 1)
    fileInput(
        inputId="controlsFCS", 
        label="Upload single-stains (FCS)", 
        accept=".fcs")
})

# render fileInput for multiplexed data
# once a spillover matrix (CSV) has been uploaded
# or single-stains have been checked
output$uploadMp <- renderUI({
    if (!is.null(vals$sm) || !is.null(vals$ffControls_metsChecked))
        tagList(
            hr(style="border-color:black"),
            fileInput(
                inputId="fcsComp", 
                label="Upload multiplexed data (FCS)", 
                multiple=TRUE,
                accept=".fcs"))
})

# toggle checkboxes
observe({
    req(input$checkbox_estSm == 1)
    updateCheckboxInput(session, "checkbox_upldSm", value=FALSE)
})
observe({
    req(input$checkbox_upldSm == 1)
    updateCheckboxInput(session, "checkbox_estSm", value=FALSE)
})

# ------------------------------------------------------------------------------
# checkboxInput: "Estimate spill from controls"
# ------------------------------------------------------------------------------

# get flowFrame of single-stained controls
ffControls <- reactive({
    ff <- vals$ffControls_metsChecked
    if (!is.null(ff)) return(ff)
    ff <- vals$ffControls_msChecked
    if (!is.null(ff)) return(ff)
    req(input$controlsFCS)
    flowCore::read.FCS(
        filename=input$controlsFCS$datapath,
        transformation=FALSE,
        truncate_max_range=FALSE)
})

# render selectInput: "Select single-positive channels"
output$selectSinglePosChs <- renderUI({
    req(input$checkbox_estSm == 1, 
        !is.null(vals$fsComp_metsChecked))
    chs <- flowCore::colnames(ffControls())
    ms <- as.numeric(CATALYST:::get_ms_from_chs(chs))
    tagList(
        hr(style="border-color:black"),
        selectSinglePosChsUI(choices=chs[!is.na(ms)]))
})

# ------------------------------------------------------------------------------
# bsButton: "Debarcode"
# ------------------------------------------------------------------------------
observe({
    test <- input$checkbox_estSm == 1 && !is.null(input$singlePosChs)
    toggleState(id="debarcodeComp", condition=test)
})

observeEvent(input$debarcodeComp, {
    disable(id="debarcodeComp")
    showNotification(h5("Assigning preliminary IDs..."), 
        duration=NULL, closeButton=FALSE, id="msg", type="message")
    ms <- as.numeric(get_ms_from_chs(input$singlePosChs))
    vals$dbFrame1Comp <- CATALYST::assignPrelim(x=ffControls(), y=ms)
    removeNotification(id="msg")
    # render deconvolution parameter UI
    output$debaParsComp <- renderUI(tagList(
        debaParsModule(module="Comp"),
        bsButton(
            inputId="compensate", 
            label="Compensate", 
            style="warning",
            size="small",
            block=TRUE)))
    enable(id="debarcodeComp")
})

# estCutoffs()
observeEvent(input$checkbox_estCutoffsComp, {
    vals$dbFrame2Comp <- CATALYST::estCutoffs(x=vals$dbFrame1Comp)
})

# toggle checkboxes
observe({
    req(input$checkbox_estCutoffsComp == 1)
    updateCheckboxInput(session, "checkbox_adjustCutoffComp", value=FALSE)
    updateCheckboxInput(session, "checkbox_globalCutoffComp", value=FALSE)
})
observe({
    req(input$checkbox_adjustCutoffComp == 1)
    updateCheckboxInput(session, "checkbox_estCutoffsComp",   value=FALSE)
    updateCheckboxInput(session, "checkbox_globalCutoffComp", value=FALSE)
})
observe({
    req(input$checkbox_globalCutoffComp == 1)
    updateCheckboxInput(session, "checkbox_estCutoffsComp",   value=FALSE)
    updateCheckboxInput(session, "checkbox_adjustCutoffComp", value=FALSE)
})

# ------------------------------------------------------------------------------
# cutoff adjustment
# ------------------------------------------------------------------------------

# toggle adjustCutoffsUI
observe({
    toggle(
        id="adjustCutoffUIComp", 
        condition=input$checkbox_adjustCutoffComp)
})

# get selectInput choices
adjustCutoffChoicesComp <- reactive({
    req(dbFrameComp())
    rownames(bc_key(dbFrameComp()))
})

# render adjustCutoffsUI
output$adjustCutoffUIComp <- renderUI({
    adjustCutoffUI(
        dbFrame=vals$dbFrame2Comp, 
        choices=adjustCutoffChoicesComp(),
        module="Comp")
})

# synchronize selectInput & numericInput w/ yield plot
observe({
    x <- selectedYieldPlotComp()
    req(!is.null(x), x != 0)
    updateSelectInput(session, "select_adjustCutoffComp", selected=x)
    updateNumericInput(session, "input_adjustCutoffComp", 
        value=paste(sep_cutoffs(vals$dbFrame2Comp)[x]))
})
observe({
    req(input$select_adjustCutoffComp)
    updateSelectInput(session, "select_yieldPlotComp", 
        selected=input$select_adjustCutoffComp)
})

# adjust cutoff upon bsButton click
observeEvent(input$button_adjustCutoffComp, {
    x <- match(input$select_adjustCutoffComp, adjustCutoffChoicesComp())
    sep_cutoffs(vals$dbFrame2Comp)[x] <- 
        as.numeric(input$input_adjustCutoffComp)
})

# ------------------------------------------------------------------------------
# global separation cutoff
# ------------------------------------------------------------------------------

# toggle globalCutoffUI
observe({
    toggle(
        id="globalCutoffUIComp", 
        condition=input$checkbox_globalCutoffComp)
})

# render globalCutoffUI
output$globalCutoffUIComp <- renderUI(globalCutoffUI(module="Comp"))

# set global cutoff upon bsButton click
observeEvent(input$button_globalCutoffComp, {
    sep_cutoffs(vals$dbFrame2Comp) <- as.numeric(input$input_globalCutoffComp)
})

# sliderInput: "Mahalanobis distance threshold"
observeEvent(input$button_mhlCutoffComp, {
    vals$mhlCutoffComp <- input$mhlCutoffComp
})

# apply cutoffs if deconvolution parameters change
dbFrameComp <- reactive({
    req(vals$dbFrame2Comp)
    applyCutoffs(vals$dbFrame2Comp, vals$mhlCutoffComp)
})

# ------------------------------------------------------------------------------
# yieldPlot()
# ------------------------------------------------------------------------------

# keep track of currently selected sample
selectedYieldPlotComp <- reactive({input$select_yieldPlotComp})

# get inputSelect choices
yieldPlotChoicesComp <- reactive({
    req(vals$dbFrame1Comp)
    ids <- rownames(bc_key(vals$dbFrame1Comp))
    setNames(c(0, ids), c("All", ids))
})

# render UI
output$yieldPlotPanelComp <- renderUI({
    req(yieldPlotChoicesComp())
    yieldPlotModule(yieldPlotChoicesComp(), "Comp")
})

# render plot
output$yieldPlotComp <- renderPlotly({
    req(dbFrameComp())
    CATALYST::plotYields(
        x=dbFrameComp(), 
        which=selectedYieldPlotComp())[[1]]
})

# renderDataTable: IDs | Counts | Cutoffs | Yields
output$summaryTblComp <- DT::renderDataTable({ 
    req(dbFrameComp())
    summary_tbl(dbFrameComp()) 
})

# next / previous buttons
observe({
    choices <- yieldPlotChoicesComp()
    n <- length(choices)
    selected <- match(selectedYieldPlotComp(), choices)
    toggleState(id="prev_yieldPlotComp", condition=selected != 1)
    toggleState(id="next_yieldPlotComp", condition=selected != n)
})
observeEvent(input$prev_yieldPlotComp, { 
    choices <- yieldPlotChoicesComp()
    selected <- match(selectedYieldPlotComp(), choices)
    updateSelectInput(session, 
        inputId="select_yieldPlotComp", 
        selected=choices[selected-1]) 
})
observeEvent(input$next_yieldPlotComp, { 
    choices <- yieldPlotChoicesComp()
    selected <- match(selectedYieldPlotComp(), choices)
    updateSelectInput(session, 
        inputId="select_yieldPlotComp", 
        selected=choices[selected+1]) 
})

# ------------------------------------------------------------------------------
# get spillover matrix, plotSpillmat()
# ------------------------------------------------------------------------------

# estimate spillover from controls
spillMatEst <- eventReactive(input$compensate, {
    req(dbFrameComp())
    CATALYST::computeSpillmat(dbFrameComp())
})

# get spillover matrix
spillMat <- reactive({
    x <- input$checkbox_upldSm
    if (!is.null(x) && x == 1 && !is.null(input$inputSpillMat)) {
        # validity check
        isCSV <- length(grep("(.csv)$", input$inputSpillMat$name)) == 1
        if (!isCSV) {
            showNotification(
                h4(strong("Input spillover matrix should be a CSV file.")),
                duration=NULL, type="error")
            return()
        }
        read.csv(input$inputSpillMat$datapath, check.names=FALSE, row.names=1)
    } else {
        spillMatEst()
    }
})

# copy original spillover matrix for adjustments
observeEvent(spillMat(), {
    vals$sm <- spillMat()
})

# plotSpillmat()
output$plotSpillmat <- renderPlot({
    req(vals$sm)
    if (input$checkbox_upldSm) {
        ms <- CATALYST:::get_ms_from_chs(colnames(vals$sm))
        CATALYST::plotSpillmat(bc_ms=ms, SM=vals$sm)
    } else if (input$checkbox_estSm) {
        ms <- rownames(bc_key(dbFrameComp()))
        CATALYST::plotSpillmat(bc_ms=ms, SM=vals$sm)
    }
})

# ------------------------------------------------------------------------------
# compensate, render download & "Go to debarcoding" button
# ------------------------------------------------------------------------------

# compensate input flowFrame(s)
fsComped <- reactive({
    req(vals$sm)
    fsApply(fsComp(), CATALYST::compCytof, vals$sm)
})

# render download buttons and bsButton "Go to debarcoding' 
# when data has been compensated
output$compensationSidebar <- renderUI({
    req(fsComped())
    compensationSidebar
})

# bsButton: "Go to debarcoding": propagate data 
# & hide FCS fileInput from debarcoding tab
observeEvent(input$goToDeba, {
    vals$keepDataComp <- TRUE
    shinyjs::hide(id="fcsDeba")
    updateTabItems(session, inputId="tabs", selected="debarcoding")
})

# ------------------------------------------------------------------------------
# before vs. after compensation scatters
# ------------------------------------------------------------------------------

# render UI
output$panel_scatters <- renderUI({
    req(spillMat())
    panel_scatters(samples=input$fcsComp$name)
})

# keep track of currently selected sample
selectedSmplComp <- reactive({
    req(input$select_smplComp)
    match(input$select_smplComp, input$fcsComp$name)
})

# ------------------------------------------------------------------------------
# next / previous sample buttons
# ------------------------------------------------------------------------------
observe({
    n <- length(fsComp())
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
    chs <- flowCore::colnames(fsComped()[[selected]])
    # default to currently selected masses 
    ms <- CATALYST:::get_ms_from_chs(chs)
    currentMs <- sapply(c(ch1(), ch2()), 
        function(i) CATALYST:::get_ms_from_chs(i))
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

# ------------------------------------------------------------------------------

# bsButton: Swap x- and y-axis
observeEvent(input$flipAxes, priority=1, {
    ch1 <- ch1()
    updateSelectInput(session, inputId="scatterCh1", selected=ch2())
    updateSelectInput(session, inputId="scatterCh2", selected=ch1)
})

# textOutput: Spillover of current interaction
output$text_spill <- renderText({
    req(ch1() %in% rownames(vals$sm) && ch2() %in% colnames(vals$sm))
    spill <- vals$sm[ch1(), ch2()]
    paste0(sprintf("%.3f", 100*spill), "%")
})

# bsButton: Adjust spill of current interaction
observe({
    toggleState(
        id="adjustSpill", 
        condition=is.numeric(input$newSpill))
})
observeEvent(input$adjustSpill, {
    vals$sm[ch1(), ch2()] <- input$newSpill/100
})

# bsButton revert current / all spill adjustments
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
    req(!is.numeric(input$cfComp), input$cfComp <= 0)
    updateNumericInput(session, inputId="cfComp", value=5)
})

observeEvent(c(input$viewCompScatter, input$flipAxes), {
    x <- ch1(); y <- ch2()
    req(x, y)
    selected <- selectedSmplComp()
    uncomped <- fsComp()[[selected]]
    comped <- fsComped()[[selected]]
    output$compScatter1 <- renderPlot(CATALYST:::plotScatter(
        es=flowCore::exprs(uncomped), x=x, y=y, cf=cfComp()))
    output$text_info1 <- renderText(text_info(uncomped, 
        isolate(input$cfComp), input$rect1, x, y))
    output$compScatter2 <- renderPlot(CATALYST:::plotScatter(
        es=flowCore::exprs(comped), x=x, y=y, cf=cfComp()))
    output$text_info2 <- renderText(text_info(comped, 
        isolate(input$cfComp), input$rect2, x, y))
})

# ------------------------------------------------------------------------------
# download handlers
# ------------------------------------------------------------------------------

# get output file names
smplNmsComp <- reactive({
    if (vals$keepDataNorm) {
        gsub(".fcs", "_comped.fcs", identifier(fsComped()))
    } else {
        req(input$fcsComp)
        gsub(".fcs", "_comped.fcs", input$fcsComp$name, ignore.case=TRUE)
    }
})

output$dwnld_comped <- downloadHandler(
    filename=function() { 
        paste0(format(Sys.Date(), "%y%m%d"), "_compensation.zip")
    },
    content=function(file) { 
        tmpdir <- tempdir()
        setwd(tmpdir)
        comped <- fsComped()
        fileNms <- smplNmsComp()
        if (input$box_setToZero) {
            lapply(seq_along(comped), function(i) {
                flowCore::exprs(comped[[i]])[
                    flowCore::exprs(comped[[i]]) < 0] <- 0
                suppressWarnings(flowCore::write.FCS(comped[[i]], fileNms[i])) 
            })
        } else {
            lapply(seq_along(comped), function(i) {
                suppressWarnings(flowCore::write.FCS(comped[[i]], fileNms[i]))
            })
        }
        zip(zipfile=file, files=fileNms)
    },
    contentType="application/zip"
)

output$dwnld_spillMat <- downloadHandler(
    filename=function() { 
        paste0(format(Sys.Date(), "%y%m%d"), "_spillMat.csv") 
    },
    content=function(file) {  
        write.csv(vals$sm, file)
    }
)