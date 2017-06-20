# read input FCS
output$debarcodingSidebar1 <- renderUI({
    if (is.null(input$fcsDeba)) return()
    vals$ffsDeba <- flowCore::read.FCS(input$fcsDeba$datapath)
    debarcodingSidebar1
})

# toggle checkboxes
observe({
    x <- input$box_csv
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_bcChs", value=FALSE)
})

observe({
    x <- input$box_bcChs
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_csv", value=FALSE)
})

# checkboxInput "Upload barcoding scheme (CSV)" & 
#               "Select single-positive channels"
output$file_csv <- renderUI({
    x <- input$box_csv
    if (is.null(x) || x == 0) return()
    fileInput(inputId="csv", label=NULL, accept=".csv")
})

output$select_bcChs <- renderUI({
    x <- input$box_bcChs
    if (is.null(x) || x == 0) return()
    selectInput(inputId="input_bcChs", label=NULL, 
                choices=flowCore::colnames(vals$ffsDeba), 
                multiple=TRUE, selectize=FALSE, size=12)
})

# actionButton: "Debarcode"
observeEvent(input$button_assignPrelim, {
    if (is.null(input$csv) & is.null(input$input_bcChs)) {
        showNotification("Please upload a barcoding scheme or select 
            single-positive channels.", type="error", closeButton=FALSE)
        return()
    }
    
    # get barcode key
    if (!is.null(input$input_bcChs)) {
        vals$DebaKey <- gsub("[[:alpha:][:punct:]]", "", input$input_bcChs)
    } else if (!is.null(input$csv)) {
        vals$DebaKey <- read.csv(input$csv$datapath, 
            check.names=FALSE, row.names=1)
    }
    
    # assign IDs
    showNotification(h5("Assigning preliminary IDs..."), id="msg", 
        type="message", duration=NULL, closeButton=FALSE)
    vals$reAssignPrelim <- CATALYST::assignPrelim(x=vals$ffsDeba, y=vals$DebaKey)
    vals$reEstCutoffs <- CATALYST::estCutoffs(x=vals$reAssignPrelim)
    removeNotification(id="msg")
    
    # inputSelect choices for yield plot and cutoff adjustment
    vals$adj_choices <- rownames(bc_key(vals$reAssignPrelim))
    vals$yp_choices  <- c(0, rownames(bc_key(vals$reAssignPrelim)))
    names(vals$yp_choices) <- c("All", rownames(bc_key(vals$reAssignPrelim)))
    
    # render yield, event and mahal plot panel
    output$yp_panel <- renderUI(yp_panel(vals$yp_choices))
    output$ep_panel <- renderUI(ep_panel(vals$ep_choices))
    output$mhl_panel <- renderUI(mhl_panel(vals$mhl_choices))
    
    output$debarcodingSidebar2 <- renderUI ({ debarcodingSidebar2 })
})

# toggle checkboxes
observe({
    x <- input$box_estCutoffs
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_adjustCutoff", value=FALSE)
    updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
    vals$reEstCutoffs <- CATALYST::estCutoffs(x=vals$reAssignPrelim)
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

# checkboxInput "Adjust population-specific cutoffs"
observe({
    x <- input$box_adjustCutoff
    toggle(id="select_adjustCutoff", condition=x)
    toggle(id="input_adjustCutoff",  condition=x)
    toggle(id="button_adjustCutoff", condition=x)
})

output$adjustCutoffs <- renderUI({
    if (input$box_adjustCutoffs == 0)
        return()
    tagList({
        div(style="display:inline-block; vertical-align:top; width:40%",
            uiOutput("select_adjustCutoff"))
        div(style="display:inline-block; width:40%",
            uiOutput("input_adjustCutoff"))
        div(style="display:inline-block; vertical-align:top",
            uiOutput("button_adjustCutoff"))
    })
})

output$select_adjustCutoff <- renderUI(
    selectInput(inputId="select_adjustCutoff", 
                label=NULL, 
                choices=vals$adj_choices))

output$input_adjustCutoff <- renderUI({
    if (input$box_adjustCutoff == 0)
        return()
    selected <- vals$adj_choices == input$select_adjustCutoff
    numericInput(inputId="input_adjustCutoff", 
                 label=NULL,
                 value=sep_cutoffs(vals$reEstCutoffs)[selected], 
                 min=0, max=1, step=.01)
})

output$button_adjustCutoff <- renderUI(
    tagList(shinyBS::bsButton(inputId="button_adjustCutoff", 
                              label=NULL,
                              icon=icon("share"),
                              style="warning"),
    shinyBS::bsTooltip(id="button_adjustCutoff",
                       title="Adjust",
                       placement="right"))
)

# synchronize selectInput w/ yield plot
observe({
    x <- input$yp_which
    if (is.null(x) || x == 0) return()
    updateSelectInput(session, "select_adjustCutoff", selected=x)
})

observeEvent(input$button_adjustCutoff, {
    selected <- vals$adj_choices == input$select_adjustCutoff
    sep_cutoffs(vals$reEstCutoffs)[selected] <- input$input_adjustCutoff
})

# checkboxInput "Enter global separation cutoff"
observe({
    x <- input$box_globalCutoff
    toggle(id="input_globalCutoff",  condition=x)
    toggle(id="button_globalCutoff", condition=x)
})

output$input_globalCutoff <- renderUI(
    numericInput(
        inputId="input_globalCutoff", 
        label=NULL, value=NULL, 
        min=0, max=1, step=.01)
)

output$button_globalCutoff <- renderUI(
    tagList(
        shinyBS::bsButton(
            inputId="button_globalCutoff", 
            label=NULL, 
            icon=icon("share"), 
            style="warning"),
        shinyBS::bsTooltip(
            id="button_globalCutoff",
            title="Apply", 
            placement="right")
    )
)

observeEvent(input$button_globalCutoff, {
    x <- input$input_globalCutoff
    if (is.null(x)) return()
    sep_cutoffs(vals$reEstCutoffs) <- as.numeric(x)
})

# sliderInput: Mahalanobis distance threshold
observeEvent(input$button_mhlCutoff, { 
    vals$mhlCutoff <- input$slider_mhlCutoff 
})

# apply cutoffs if deconvolution parameters change
observe({
    x <- vals$reEstCutoffs
    if (is.null(x)) return()
    vals$reApplyCutoffs <- CATALYST::applyCutoffs(x, vals$mhlCutoff)
})

# ------------------------------------------------------------------------------
# yield, event and mahal plot
# ------------------------------------------------------------------------------

observe({
    x <- vals$reApplyCutoffs
    if (is.null(x)) return()
    output$plotYields <- renderPlot(
        plotYields(x, input$yp_which))
    output$plotEvents <- renderPlot(
        plotEvents(x, input$ep_which, as.numeric(input$n_events)))
    output$plotMahal <- renderPlot(
        plotMahal(x, input$mhl_which, input$input_mhlCofactor))
})

# inputSelect choices for event and mahal plot
observe({
    if (is.null(vals$reApplyCutoffs)) return()
    vals$ids <- sort(unique(bc_ids(vals$reApplyCutoffs)))
    vals$ep_choices  <- vals$ids
    vals$mhl_choices <- vals$ids[vals$ids != 0]
})

# summary table: IDs | Counts | Cutoffs | Yields
output$table_summary <- DT::renderDataTable({ 
    if (is.null(vals$reApplyCutoffs)) return()
    summary_tbl(vals$reApplyCutoffs) 
})

# ------------------------------------------------------------------------------
# next / previous buttons
# ------------------------------------------------------------------------------
observeEvent(input$yp_prev, { 
    x <- which(vals$yp_choices == input$yp_which)
    if (x == 1) return()
    updateSelectInput(session, "yp_which", selected=vals$yp_choices[x-1]) 
})

observeEvent(input$yp_next, { 
    x <- which(vals$yp_choices == input$yp_which)
    if (x == length(vals$yp_choices)) return()
    updateSelectInput(session, "yp_which", selected=vals$yp_choices[x+1]) 
})

observeEvent(input$ep_prev, { 
    x <- which(vals$ep_choices == input$ep_which)
    if (x == 1) return()
    updateSelectInput(session, "ep_which", selected=vals$ep_choices[x-1]) 
})

observeEvent(input$ep_next, { 
    x <- which(vals$ep_choices == input$ep_which)
    if (x == length(vals$ep_choices)) return()
    updateSelectInput(session, "ep_which", selected=vals$ep_choices[x+1])
})

observeEvent(input$mhl_prev, { 
    x <- which(vals$mhl_choices == input$mhl_which)
    if (x == 1) return()
    updateSelectInput(session, "mhl_which", selected=vals$mhl_choices[x-1])
})

observeEvent(input$mhl_next, { 
    x <- which(vals$mhl_choices == input$mhl_which)
    if (x == length(vals$mhl_choices)) return()
    updateSelectInput(session, "mhl_which", selected=vals$mhl_choices[x+1])
})

# ------------------------------------------------------------------------------
# synchronize plots
# ------------------------------------------------------------------------------
observe({
    x <- input$yp_which
    if (is.null(x) || x == 0) return() 
    updateSelectInput(session, "ep_which",  selected=x)
    updateSelectInput(session, "mhl_which", selected=x)
})

observe({
    x <- input$ep_which
    if (is.null(x) || x == 0) return() 
    updateSelectInput(session, "yp_which",  selected=x)
    updateSelectInput(session, "mhl_which", selected=x)
})

observe({
    x <- input$mhl_which
    if (is.null(x)) return()
    updateSelectInput(session, "yp_which",  selected=x)
    updateSelectInput(session, "ep_which",  selected=x)
})
