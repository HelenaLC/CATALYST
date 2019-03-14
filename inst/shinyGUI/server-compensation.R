# ==============================================================================
# COMPENSATION
# ==============================================================================

# get flowSet to compensated: 
# initially, use normalized samples or read input FCS file(s)
# once duplicate masses or metal mismatches where checked,
# the newly constructed data will be used
fsComp <- reactive({
    if (!is.null(fs <- vals$fsComp_metsChecked)) {
        return(fs)
    } else if (!is.null(fs <- vals$fsComp_msChecked)) {
        # if provided, add controls to flowSet
        if (!is.null(ffControls())) {
            controls <- ffControls()
            description(controls)[c("GUID", "ORIGINALGUID")] <- 
                identifier(controls) <- input$fcsControls$name
            fs <- as(c(controls, lapply(seq_along(fs), 
                function(i) fs[[i]])), "flowSet")
        }
        return(fs)
    } else if (vals$keepDataNorm) {
            ffs <- ffsNormed()
            if (input$box_removeBeads) {
                for (i in seq_along(ffs)) {
                    removed <- mhlDists()[[i]] < vals$mhlCutoffs[i]
                    ffs[[i]] <- ffs[[i]][!removed, ]
                }
            }
            fs <- as(ffs, "flowSet")
            nms <- gsub(".fcs", "_normed.fcs", smplNmsNorm(), ignore.case=TRUE)
    } else {
        # read & check input FCS files
        fcs <- req(input$fcsComp)
        n <- nrow(fcs)
        req(check_FCS_fileInput(fcs, n))
        ffs <- lapply(seq_len(n), function(i) 
            flowCore::read.FCS(
                filename=fcs[[i, "datapath"]],
                transformation=FALSE,
                truncate_max_range=FALSE))
        fs <- as(ffs, "flowSet")
        nms <- input$fcsComp$name
    }
    for (i in seq_along(fs))
        description(fs[[i]])[c("GUID", "ORIGINALGUID")] <- 
            identifier(fs[[i]]) <- nms[i]
    return(fs)
})

# render fileInput for spillover matrix CSV or single-stains FCS
# according to radioButtons selection
output$upload_or_est_sm_UI <- renderUI({
    switch(input$upload_or_est_sm, 
        upload_sm=fileInput(
            inputId="input_sm", 
            label="Upload spillover matrix (CSV)", 
            accept=".csv"),
        est_sm=fileInput(
            inputId="fcsControls", 
            label="Upload single-stains (FCS)", 
            accept=".fcs"))
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
    req(input$fcsControls)
    # check validity of input FCS file
    valid <- check_FCS_fileInput(input$fcsControls, 1)
    if (!valid) return()
    flowCore::read.FCS(
        filename=input$fcsControls$datapath,
        transformation=FALSE,
        truncate_max_range=FALSE)
})

# render selectInput: "Select single-positive channels"
output$selectSinglePosChs <- renderUI({
    req(input$upload_or_est_sm == "est_sm", 
        !is.null(vals$fsComp_metsChecked))
    chs <- flowCore::colnames(ffControls())
    ms <- as.numeric(CATALYST:::.get_ms_from_chs(chs))
    tagList(
        hr(style="border-color:black"),
        selectSinglePosChsUI(choices=chs[!is.na(ms)]))
})

# ------------------------------------------------------------------------------
# bsButton: "Debarcode"
# ------------------------------------------------------------------------------
observe({
    test <- input$upload_or_est_sm == "est_sm" && !is.null(input$singlePosChs)
    toggleState(id="debarcodeComp", condition=test)
})

observeEvent(input$debarcodeComp, {
    disable(id="debarcodeComp")
    showNotification(h4(strong("Assigning preliminary IDs...")), 
        duration=NULL, closeButton=FALSE, id="msg", type="message")
    ms <- as.numeric(CATALYST:::.get_ms_from_chs(input$singlePosChs))
    vals$dbFrame1Comp <- CATALYST::assignPrelim(x=ffControls(), y=ms)
    removeNotification(id="msg")
    # estCutoffs()
    showNotification("Estimating separation cutoffs...",
        id="estimating_sep_cutoffs", duration=NULL, closeButton=FALSE)
    vals$cutoff_ests_comp <- sep_cutoffs(CATALYST::estCutoffs(x=vals$dbFrame1Comp))
    removeNotification(id="estimating_sep_cutoffs")
    # render deconvolution parameter UI
    output$debaParsComp <- renderUI(
        debaParsModule(module="Comp"))
    output$est_sm <- renderUI(
        bsButton(
            inputId="est_sm", 
            label="Estimate spillover", 
            style="warning",
            size="small",
            block=TRUE))
    enable(id="debarcodeComp")
})

# ------------------------------------------------------------------------------
# cutoff adjustment
# ------------------------------------------------------------------------------
# render UI for cutoff adjustment or input for global cutoff
output$deba_cutoffs_UIComp <- renderUI({
    switch(input$deba_cutoffsComp,
        est_cutoffs=NULL,
        adj_cutoffs=
            adjustCutoffUI(
                dbFrame=vals$dbFrame2Comp, 
                choices=adjustCutoffChoicesComp(),
                module="Comp"),
        global_cutoff=
            globalCutoffUI(module="Comp"))
})
# use cutoff estimates
observeEvent(input$deba_cutoffsComp == "est_cutoffs", ignoreInit=TRUE, {
    sep_cutoffs(vals$dbFrame1Comp) <- vals$cutoff_ests_comp
    vals$dbFrame2Comp <- vals$dbFrame1Comp
})

# get selectInput choices for cutoff adjustment
adjustCutoffChoicesComp <- reactive({
    req(dbFrameComp())
    rownames(bc_key(dbFrameComp()))
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

# set global cutoff upon bsButton click
observeEvent(input$button_globalCutoffComp, {
    sep_cutoffs(vals$dbFrame2Comp) <- input$input_globalCutoffComp
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
estimated_sm <- eventReactive(input$est_sm, {
    req(dbFrameComp())
    showNotification(h4(strong("Estimating spillover...")),
        id="msg", type="message", duration=NULL, closeButton=FALSE)
    sm <- CATALYST::computeSpillmat(dbFrameComp())
    removeNotification(id="msg")
    return(sm)
})

input_sm <- eventReactive(input$input_sm, {
    is_csv <- length(grep("(.csv)$", input$input_sm$name)) == 1
    if (!is_csv) {
        showNotification(type="error", closeButton=FALSE,
            h4(strong("Input spillover matrix should be a CSV file.")))
        return()
    }
    sm <- read.csv(input$input_sm$datapath, row.names=1)
    sm <- tryCatch(CATALYST:::.check_sm(sm, vals$isotope_list), error=function(e) e)
    if (inherits(sm, "error")) {
        showNotification(type="error", closeButton=FALSE,
            h4(strong("Input spillover matrix seems to be invalid.")))
        return()
    } 
    return(sm)
})

# get spillover matrix
sm <- reactive({
    x <- input$upload_or_est_sm
    switch(x,
        upload_sm = req(input_sm()),
        est_sm = req(estimated_sm()))
})

# copy original spillover matrix for adjustments
observeEvent(sm(), { 
    vals$sm <- sm()
})

# plotSpillmat()
output$plotSpillmat <- renderPlotly({
    req(vals$sm)
    x <- input$upload_or_est_sm
    ms <- switch (x,
        upload_sm = CATALYST:::.get_ms_from_chs(colnames(vals$sm)),
        est_sm = rownames(bc_key(dbFrameComp())))
    pdf(NULL)
    CATALYST::plotSpillmat(bc_ms=ms, SM=vals$sm, isotope_list=vals$isotope_list)
})

# ------------------------------------------------------------------------------
# compensate, render download & "Go to debarcoding" button
# ------------------------------------------------------------------------------

# toggle bsButton "Compensate"
observe({
    test <- req(sm(), fsComp())
    toggleState(id="compensate", condition=test)
})

# compensate input flowFrame(s)
nnlsComped <- reactive({
    req(fsComp(), vals$sm, input$compensation_method == "nnls")
    showNotification(h4(strong("Compensating using NNLS...")), 
        id="msg", type="message", duration=NULL, closeButton=FALSE)
    fs <- CATALYST::compCytof(fsComp(), vals$sm, NULL, "nnls", vals$isotope_list)
    removeNotification(id="msg")
    return(fs)
})
flowComped <- reactive({
    req(input$compensation_method == "flow", vals$sm, fsComp())
    showNotification(h4(strong("Compensating using method 'flow'...")), 
        id="msg", type="message", duration=NULL, closeButton=FALSE)
    fs <- CATALYST::compCytof(fsComp(), vals$sm, NULL, "flow", vals$isotope_list)
    removeNotification(id="msg")
    return(fs)
})
fsComped <- eventReactive(input$compensate, {
    disable(id="compensate")
    method <- input$compensation_method
    fs <- switch(method, nnls=nnlsComped(), flow=flowComped())
    nms <- keyword(fs, "ORIGINALGUID")
    for (i in seq_along(fs))
        description(fs[[i]])$GUID <- identifier(fs[[i]]) <- nms[i]
    enable(id="compensate")
    return(fs)
})

# render download buttons 
# & bsButton "Go to debarcoding" 
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
observeEvent(input$compensate, {
    output$panel_scatters <- renderUI(
        panel_scatters(samples=smplNmsComp()))
})
#  <- renderUI({
#     req(sm())
#     panel_scatters(samples=smplNmsComp())
# })

# keep track of currently selected sample
selectedSmplComp <- reactive({
    req(input$select_smplComp)
    match(input$select_smplComp, smplNmsComp())
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
        inputId="select_smplComp", 
        selected=smplNmsComp()[selectedSmplComp()-1]) 
})
observeEvent(input$next_smplComp, { 
    updateSelectInput(session, 
        inputId="select_smplComp", 
        selected=smplNmsComp()[selectedSmplComp()+1]) 
})

# update channel selection choices upon sample change
observe({
    selected <- selectedSmplComp()
    req(selected)
    chs <- flowCore::colnames(fsComped()[[selected]])
    # default to currently selected masses 
    ms <- CATALYST:::.get_ms_from_chs(chs)
    currentMs <- sapply(c(ch1(), ch2()), 
        function(i) CATALYST:::.get_ms_from_chs(i))
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
    vals$sm[ch1(), ch2()] <- sm()[ch1(), ch2()]
})
observeEvent(input$revertAll, {
    vals$sm <- sm()
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

# downsample to 5e3 events per sample
inds <- eventReactive(fsComp(), {
    nEvents <- fsApply(fsComp(), nrow)
    n <- pmin(nEvents, 5e3)
    lapply(seq_along(n), function(i) sample(nEvents[i], n[i]))
})

# subset raw & compensated data
es0 <- reactive({
    req(inds())
    lapply(seq_along(fsComp()), function(i)
        exprs(fsComp()[[i]][inds()[[i]], ]))
})
esC <- reactive({
    req(fsComped(), inds())
    lapply(seq_along(fsComped()), function(i)
        exprs(fsComped()[[i]][inds()[[i]], ]))
})

# render scatter plots and text outputs
output$compScatter1 <- renderPlot({
    req(es0(), ch1(), ch2(), cfComp())
    exprs <- es0()[[selectedSmplComp()]]
    CATALYST:::.plotScatter(exprs, ch1(), ch2(), cfComp())
})
output$compScatter2 <- renderPlot({
    req(esC(), ch1(), ch2(), cfComp())
    exprs <- esC()[[selectedSmplComp()]]
    CATALYST:::.plotScatter(exprs, ch1(), ch2(), cfComp())
})
output$text_info1 <- renderText({
    req(fsComp(), selectedSmplComp(), 
        cfComp(), input$rect1, ch1(), ch2())
    ff <- fsComp()[[selectedSmplComp()]]
    text_info(ff, cfComp(), input$rect1, ch1(), ch2())
})
output$text_info2 <- renderText({
    req(fsComped(), selectedSmplComp(), 
        cfComp(), input$rect1, ch1(), ch2())
    ff <- fsComped()[[selectedSmplComp()]]
    text_info(ff, cfComp(), input$rect2, ch1(), ch2())
})

# ------------------------------------------------------------------------------
# download handlers
# ------------------------------------------------------------------------------

# get sample and output file names
smplNmsComp <- reactive({
    req(fsComp())
    keyword(fsComp(), "ORIGINALGUID")
})
outNmsComp <- reactive({
    req(smplNmsComp())
    gsub(".fcs", "_comped.fcs", smplNmsComp(), ignore.case=TRUE)
})

output$dwnld_comped <- downloadHandler(
    filename=function() { 
        paste0(format(Sys.Date(), "%y%m%d"), "_compensation.zip")
    },
    content=function(file) { 
        tmpdir <- tempdir()
        setwd(tmpdir)
        comped <- fsComped()
        fileNms <- outNmsComp()
        lapply(seq_along(comped), function(i) {
            suppressWarnings(flowCore::write.FCS(comped[[i]], fileNms[i]))
        })
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
