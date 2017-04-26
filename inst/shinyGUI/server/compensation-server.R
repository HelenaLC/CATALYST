# read input FCS file(s)
observe({
    x <- input$fcsComp
    if (is.null(x)) return()
    for (i in seq_len(nrow(x)))
        vals$ffsComp[[i]] <- flowCore::read.FCS(x[[i, "datapath"]])
    output$compensationSidebar1 <- renderUI({ compensationSidebar1 })
})

observeEvent(input$fcsComp, {
    if (!is.null(input$input_bcChs) & !is.null(vals$re2)) {
        updateCheckboxInput(session,
                            inputId="box_estSm",
                            value=TRUE)
        # auto-estimate spillover matrix if possible
        vals$sm <- CATALYST::computeSpillmat(x=vals$re2)
    } else {
        updateCheckboxInput(session,
                            inputId="box_upldSm",
                            value=TRUE)
    }
})

# ------------------------------------------------------------------------------
# toggle checkboxes
# ------------------------------------------------------------------------------

observe({
    x <- input$box_upldSm
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_estSm", value=FALSE)
})

observe({
    x <- input$box_estSm
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_upldSm", value=FALSE)
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
# checkboxInput "Estimate spill from single-stained controls"
# ------------------------------------------------------------------------------

output$enterTrim <- renderUI({
    if (input$box_estSm == 1 & !is.null(input$input_bcChs) & !is.null(vals$re2))
        enterTrim
})

output$panel_estTrim <- renderUI({
    x <- input$box_estSm
    if (!is.null(x) & x == 1 & !is.null(input$input_bcChs) & !is.null(vals$re2))
        panel_estTrim
})

observe({
    toggleState(id="button_estTrim", 
                condition=(is.numeric(input$estTrim_min) 
                           & is.numeric(input$estTrim_max) 
                           & is.numeric(input$estTrim_stp) 
                           & input$estTrim_stp != 0))
    toggleState(id="button_enterTrim", 
                condition=is.numeric(input$input_enterTrim))
})

# actionButton "Go"
observeEvent(input$button_estTrim, {
    min <- input$estTrim_min
    max <- input$estTrim_max
    stp <- input$estTrim_stp
    l <- length(seq(min, max, stp))
    showNotification(paste("Estimating spill for", l, "trim values..."),
                     id="msg", type="message", duration=NULL, closeButton=FALSE)
    p <- CATALYST::estTrim(vals$re2, min, max, stp)
    output$plot_estTrim <- renderPlotly(p)
    removeNotification(id="msg")
})

# actionButton "Estimate spillover"
observeEvent(input$button_enterTrim, {
    showNotification(h5("Estimating spillover..."), id="msg", 
                     type="message", duration=NULL, closeButton=FALSE)
    vals$trm <- input$input_enterTrim
    vals$sm  <- CATALYST::computeSpillmat(x=vals$re2, trim=vals$trm)
    removeNotification(id="msg")
})

# ------------------------------------------------------------------------------
# checkboxInput "Upload spillover matrix (CSV)"
# ------------------------------------------------------------------------------

output$inputSm <- renderUI({
    if (input$box_upldSm)
        fileInput("inputSm", NULL, accept=".csv")
})

observe({
    if (!is.null(input$inputSm))
        vals$sm0 <- vals$sm <- as.matrix(
            read.csv(input$inputSm$datapath, 
                     check.names=FALSE, row.names=1))
})

# ------------------------------------------------------------------------------
# render panels for spillover matrix heatmap 
# & before vs. after compensation scatters
# when spillover matrix is available
observeEvent(vals$sm, once=TRUE, {
    output$panel_plotSpillmat <- renderUI(panel_plotSpillmat)
    # default to interaction with highest spill
    sm <- isolate(vals$sm)
    inds <- which(sm == max(sm[sm != 1]), TRUE)
    chs <- flowCore::colnames(vals$ffsComp[[1]])
    ms <- gsub("[[:punct:][:alpha:]]", "", chs)
    msSm <- lapply(dimnames(sm), function(i) 
        gsub("[[:punct:][:alpha:]]", "", i))
    x <- match(msSm[[1]][inds[1]], ms)
    y <- match(msSm[[2]][inds[2]], ms)
    output$panel_scatters <- renderUI(
        panel_scatters(samples=input$fcsComp$name,
                       channels=chs, x=x, y=y))        
})

# keep track of currently selected sample
selectedComp <- reactive({
    if (!is.null(input$select_sampleComp)) {
        which(input$fcsComp$name == input$select_sampleComp)
    } else {
        NULL
    }
})

# next / previous sample buttons
observe({
    n <- length(vals$ffsComp)
    toggleState(id="prevSmplComp", condition=selectedComp()!=1)
    toggleState(id="nextSmplComp", condition=selectedComp()!=n)
})

observeEvent(input$prevSmplComp, { 
    updateSelectInput(session, "select_sampleComp", 
                      selected=input$fcsComp$name[selectedComp()-1]) 
})

observeEvent(input$nextSmplComp, { 
    updateSelectInput(session, "select_sampleComp", 
                      selected=input$fcsComp$name[selectedComp()+1]) 
})

# update channel selection choices upon sample change
observe({
    if (is.null(selectedComp())) return()
    chs <- flowCore::colnames(isolate(vals$ffsComped)[[selectedComp()]])
    # default to currently selected masses 
    ms <- gsub("[[:punct:][:alpha:]]", "", chs)
    currentMs <- sapply(c(isolate(input$scatterCh1), isolate(input$scatterCh2)), 
                        function(i) gsub("[[:punct:][:alpha:]]", "", i))
    inds <- match(currentMs, ms)
    # if non-existent, plot first 2 mass channels
    if (any(is.na(inds))) 
        inds <- which(!is.na(as.numeric(ms)))[1:2]
    updateSelectInput(session, "scatterCh1", choices=chs, selected=chs[inds[1]])
    updateSelectInput(session, "scatterCh2", choices=chs, selected=chs[inds[2]])
})

# actionButton: Swap x- and y-axis
observeEvent(input$flipAxes, {
    ch1 <- input$scatterCh1
    updateSelectInput(session,
                      inputId="scatterCh1",
                      selected=input$scatterCh2)
    updateSelectInput(session,
                      inputId="scatterCh2",
                      selected=ch1)
})

# textOutput: Spillover of current interaction
output$text_spill <- renderText({
    #(input$flipAxes)
    paste0(sprintf("%.3f", 100*vals$sm[
        isolate(input$scatterCh1), input$scatterCh2]), "%")
})

# actionButton: Adjust spill of current interaction
observe(toggleState(id="adjustSpill", condition=is.numeric(input$newSpill)))
observeEvent(input$adjustSpill, {
    vals$sm[input$scatterCh1, input$scatterCh2] <- input$newSpill/100
})

# actionButtons revert current / all spill adjustments
observeEvent(input$revert, {
    vals$sm[input$scatterCh1, input$scatterCh2] <-
        vals$sm0[input$scatterCh1, input$scatterCh2]
})

observeEvent(input$revertAll, {
    vals$sm <- vals$sm0
})

# render spillover matrix heatmap 
output$plot_plotSpillmat <- renderPlot({
    if (is.null(vals$sm)) return()
    input$button_newSpill
    if (input$box_estSm == 1) {
        CATALYST::plotSpillmat(bc_ms=vals$key, SM=isolate(vals$sm))
    } else {
        ms <- gsub("[[:punct:][:alpha:]]", "", colnames(vals$sm))
        CATALYST::plotSpillmat(bc_ms=ms[!is.na(as.numeric(ms))], 
                               SM=isolate(vals$sm))
    }
})

# compensate input flowFrame(s)
observe({
    if (is.null(vals$sm)) return()
    vals$ffsComped <- lapply(vals$ffsComp, CATALYST::compCytof, vals$sm)
    if (!is.null(vals$ff1) & !is.null(input$input_bcChs)) {
        output$text_compCytof_1 <- renderText(
            "WARNING: Compensation is likely to be inaccurate.
            Spill values for the following interactions
            have not been estimated:")
        output$text_compCytof_2 <- renderPrint(
            cat(capture.output(type="message", 
                               vals$ffsComped <- c(compCytof(x=vals$ff1, y=vals$sm),
                                                   vals$ffsComped))[-c(1:3)], sep="\n"))
    }
})

# ------------------------------------------------------------------------------
# before vs. after compensation scatters
# ------------------------------------------------------------------------------

# default to cofactor 5 if input is invalid
observe({
    if (!is.numeric(input$cfComp) || input$cfComp == 0)
        updateNumericInput(session, inputId="cfComp", value=5)
})

# left scatter (uncompensated)
output$compScatter1 <- renderPlot({ 
    (input$button_view )#|| input$button_newSpill)
    x <- vals$ffsComp[[selectedComp()]]
    if (is.null(x) || !is.numeric(input$cfComp)) return()
    plotScatter(es=exprs(x), n=25e3,
                x=isolate(input$scatterCh1), 
                y=input$scatterCh2, 
                cf=input$cfComp)
})
output$text_info1 <- renderText({ 
    text_info(
        vals$ffsComp[[selectedComp()]], 
        isolate(input$cfComp), 
        input$rect1, 
        isolate(input$scatterCh1), 
        isolate(input$scatterCh2))
})   

# right scatter (compensated)
output$compScatter2 <- renderPlot({
    (input$button_view )#| input$button_newSpill)
    x <- vals$ffsComped[[selectedComp()]]
    if (is.null(x) || !is.numeric(input$cfComp)) return()
    plotScatter(es=exprs(x), n=25e3,
                x=isolate(input$scatterCh1), 
                y=input$scatterCh2, 
                cf=input$cfComp)
})
output$text_info2 <- renderText({ 
    input$button_view 
    text_info(
        vals$ffsComped[[selectedComp()]], 
        isolate(input$cfComp), 
        input$rect2, 
        isolate(input$scatterCh1), 
        isolate(input$scatterCh2))
})

# ------------------------------------------------------------------------------

# check that all input flowFrames have been compensated
# & render download button
output$compensationSidebar2 <- renderUI({
    if (any(vapply(vals$ffsComped, is.null, logical(1))))
        return()
    compensationSidebar2
})