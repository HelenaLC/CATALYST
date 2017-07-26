# ==============================================================================
# COMPENSATION
# ==============================================================================

# get input flowFrames: 
# use debarcoded samples or read input FCS files
ffsComp <- reactive({
    if (vals$alterPars) {
        vals$newMp
    } else if (vals$keepDataDeba) {
        ids <- sort(unique(bc_ids(dbFrameDeba())))
        ids <- ids[ids != 0]
        vals$compFiles <- TRUE
        lapply(ids, function(id) 
            ffDeba()[bc_ids(dbFrameDeba()) == id, ])
    } else {
        req(input$fcsComp) 
        vals$compFiles <- TRUE
        lapply(seq_len(nrow(input$fcsComp)), function(i) 
            flowCore::read.FCS(
                filename=input$fcsComp[[i, "datapath"]],
                transformation=FALSE,
                truncate_max_range=FALSE))
    }
})

# render FCS file editing modal
observe({
    req(chInds(), !vals$checkPars)
    inds <- chInds()
    showModal(editFCS())
    output$chList1 <- renderUI(chList(ffControls(),   inds[[1]], "ss"))
    output$chList2 <- renderUI(chList(ffsComp()[[1]], inds[[2]], "mp"))
    vals$checkPars <- TRUE
})

# get indicies of masses matching between
# multiplexed and single-stained flowFrames
chInds <- reactive({
    req(ffsComp(), ffControls())
    chs1 <- flowCore::colnames(ffControls())
    chs2 <- flowCore::colnames(ffsComp()[[1]])
    ms1 <- gsub("[[:punct:][:alpha:]]", "", chs1)
    ms2 <- gsub("[[:punct:][:alpha:]]", "", chs2)
    ex1 <- is.na(as.numeric(ms1))
    ex2 <- is.na(as.numeric(ms2))
    ms <- intersect(ms1[!ex1], ms2[!ex2])
    inds1 <- ms1 %in% ms
    inds2 <- ms2 %in% ms
    list(inds1, inds2)
})

# color textInputs depending on (mis)match between metals
observe({
    req(input$ss1)
    inds <- chInds()
    chs1 <- flowCore::colnames(ffControls())  [inds[[1]]]
    chs2 <- flowCore::colnames(ffsComp()[[1]])[inds[[2]]]
    mets1 <- gsub("[[:punct:][:digit:]+Di]", "", chs1)
    mets2 <- gsub("[[:punct:][:digit:]+Di]", "", chs2)
    diff <- mets1 != mets2
    n <- length(diff)
    cols <- rep("#cce698", n)
    cols[diff] <- "#ffcccc"
    for (i in 1:2) {
        id <- c("ss", "mp")[i] 
        for (j in seq_len(n)) {
            js$textInputCol(id=paste0(id, j), col=cols[j])
        }
    }
})

# bsButton: "Done":
# alter parameter names of input flowFrames
observeEvent(input$done, {
    inds <- chInds()
    chs1 <- flowCore::colnames(ffControls())  [inds[[1]]]
    chs2 <- flowCore::colnames(ffsComp()[[1]])[inds[[2]]]
    mets1 <- gsub("[[:punct:][:digit:]+Di]", "", chs1)
    mets2 <- gsub("[[:punct:][:digit:]+Di]", "", chs2)
    diff <- mets1 != mets2
    n <- length(diff)
    # single-stained
    ss <- ffControls()
    pars <- flowCore::colnames(ss)
    pars[inds[[1]]] <- sapply(seq_len(n), 
        function(i) input[[paste0("ss", i)]])
    vals$newSs <- new_ff(
        data=flowCore::exprs(ss),
        pars=pars,
        desc=ss@parameters@data$desc)
    # multiplexed
    mp <- ffsComp()
    vals$newMp <- 
        lapply(seq_along(mp), function(i) {
            ff <- mp[[i]]
            pars <- flowCore::colnames(ff)
            pars[inds[[2]]] <- sapply(seq_len(n), 
                function(j) input[[paste0("mp", j)]])
            new_ff(
                data=flowCore::exprs(ff),
                pars=pars,
                desc=ff@parameters@data$desc)
        })
    vals$alterPars <- TRUE
    removeModal()
})

# expand sidebar once data has been uploaded
output$compensationSidebar1 <- renderUI({
    req(isTRUE(vals$compFiles))
    compensationSidebar1
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
# spillover matrix selection & plotting
# ------------------------------------------------------------------------------

# render fileInput for single-stained controls
output$uploadControls <- renderUI({
    req(input$checkbox_estSm == 1)
    fileInput(
            inputId="controlsFCS", 
            label="Upload FCS", 
            accept=".fcs")
})

# get flowFrame of single-stained controls
ffControls <- reactive({
    if (vals$alterPars) {
        vals$newSs
    } else {
        req(input$controlsFCS)
        flowCore::read.FCS(
            filename=input$controlsFCS$datapath,
            transformation=FALSE,
            truncate_max_range=FALSE)
    }
})

# render selectInput and bsButton
output$selectSinglePosChs <- renderUI({
    req(input$checkbox_estSm == 1, ffControls())
    chs <- flowCore::colnames(ffControls())
    ms <- as.numeric(gsub("[[:alpha:][:punct:]]", "", chs))
    selectSinglePosChsUI(choices=chs[!is.na(as.numeric(ms))])
})

# ------------------------------------------------------------------------------
# bsButton: "Debarcode"
# ------------------------------------------------------------------------------
observe({
    test <- input$checkbox_estSm == 1 && !is.null(input$singlePosChs)
    toggleState(id="debarcodeComp", condition=test)
})

observeEvent(input$debarcodeComp, {
    showNotification(h5("Assigning preliminary IDs..."), 
        duration=NULL, closeButton=FALSE, id="msg", type="message")
    vals$dbFrame1Comp <- CATALYST::assignPrelim(
        x=ffControls(), 
        y=as.numeric(gsub("[[:alpha:][:punct:]]", "", input$singlePosChs)))
    removeNotification(id="msg")
    # render deconvolution parameter UI
    output$debaParsComp <- renderUI(debaParsModule(module="Comp"))
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
output$yieldPlotComp <- renderPlot({
    req(dbFrameComp())
    CATALYST::plotYields(
        x=dbFrameComp(), 
        which=selectedYieldPlotComp())
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
# get spillover matrix & plotSpillmat()
# ------------------------------------------------------------------------------

# render fileInput for spillover matrix CSV
output$uploadSpillMat <- renderUI({
    req(input$checkbox_upldSm == 1)
    fileInput(
        inputId="inputSpillMat", 
        label=NULL, 
        accept=".csv")
})

# get spillover matrix
spillMat <- reactive({
    req(ffsComp())
    x <- input$checkbox_upldSm
    if (!is.null(x) && x == 1 && !is.null(input$inputSpillMat)) {
        read.csv(input$inputSpillMat$datapath, check.names=FALSE, row.names=1)
    } else {
        x <- input$checkbox_estSm
        if (!is.null(x) && x == 1 && !is.null(dbFrameComp()))
            computeSpillmat(x=dbFrameComp())
    }
})

# store original spillover matrix
observeEvent(spillMat(), {
    vals$sm <- spillMat()
})

# plotSpillmat()
output$plotSpillmat <- renderPlot({
    req(vals$sm)
    if (input$checkbox_upldSm) {
        ms <- gsub("[[:punct:][:alpha:]]", "", colnames(vals$sm))
        CATALYST::plotSpillmat(bc_ms=ms, SM=vals$sm)
    } else if (input$checkbox_estSm) {
        ms <- rownames(bc_key(dbFrameComp()))
        CATALYST::plotSpillmat(bc_ms=ms, SM=vals$sm)
    }
})

# compensate input flowFrame(s)
ffsComped <- reactive({
    req(vals$sm)
    print(is.list(ffsComp()))
    print(compCytof(ffsComp()[[1]], vals$sm))
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

# get output file names
smplNmsComp <- reactive({
    if (vals$keepDataDeba) {
        dbFrame <- dbFrameDeba()
        ids <- rownames(bc_key(dbFrame))
        paste0(smplNmsDeba(), "_comped.fcs")[
            c(0, ids) %in% sort(unique(bc_ids(dbFrame)))]
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
        ffs <- ffsComped()
        fileNms <- smplNmsComp()
        if (input$box_setToZero) {
            lapply(seq_along(ffs), function(i) {
                flowCore::exprs(ffs[[i]])[flowCore::exprs(ffs[[i]]) < 0] <- 0
                suppressWarnings(flowCore::write.FCS(ffs[[i]], fileNms[i])) 
            })
        } else {
            lapply(seq_along(ffs), function(i) {
                suppressWarnings(flowCore::write.FCS(ffs[[i]], fileNms[i]))
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