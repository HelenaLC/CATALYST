# ------------------------------------------------------------------------------
# toggle checkboxes
# ------------------------------------------------------------------------------
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

# "Select single-positive channels"
output$select_bcChs <- renderUI({
    x <- input$box_bcChs
    if (is.null(x) || x == 0) return()
    selectInput("input_bcChs", NULL, choices=flowCore::colnames(vals$ff1), 
        multiple=TRUE, selectize=FALSE, size=12)
})

# ------------------------------------------------------------------------------
# toggle checkboxes
# ------------------------------------------------------------------------------
observe({
    x <- input$box_estCutoffs
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_adjustCutoff", value=FALSE)
    updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
    vals$re1 <- CATALYST::estCutoffs(x=vals$re0)
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
# "Adjust population-specific cutoffs"
# ------------------------------------------------------------------------------
output$select_adjustCutoff <- renderUI({
    x <- input$box_adjustCutoff
    if (is.null(x) || x == 0) return()
    selectInput("select_adjustCutoff", NULL, choices=vals$adj_choices)
})

# synchronize selectInput w/ yield plot
observe({
    x <- input$yp_which
    if (is.null(x) || x == 0) return()
    updateSelectInput(session, "select_adjustCutoff", selected=input$yp_which)
})

output$input_adjustCutoff <- renderUI({
    x <- input$box_adjustCutoff
    if (is.null(x) || x == 0) return()
    selected <- vals$adj_choices == input$select_adjustCutoff
    numericInput("input_adjustCutoff", NULL, 
        value=sep_cutoffs(vals$re1)[selected], min=0, max=1, step=.01)
})

output$button_adjustCutoff <- renderUI({
    x <- input$box_adjustCutoff
    if (is.null(x) || x == 0) return()
    actionButton("button_adjustCutoff", "Adjust", style=ylw_button)
})

observeEvent(input$button_adjustCutoff, {
    selected <- vals$adj_choices == input$select_adjustCutoff
    sep_cutoffs(vals$re1)[selected] <- input$input_adjustCutoff
})

# ------------------------------------------------------------------------------
# "Enter global separation cutoff"
# ------------------------------------------------------------------------------
output$input_globalCutoff <- renderUI({
    x <- input$box_globalCutoff
    if (is.null(x) || x == 0) return()
    numericInput("input_globalCutoff", NULL, NULL, min=0, max=1, step=.01)
})

observeEvent(input$button_globalCutoff, {
    x <- input$input_globalCutoff
    if (is.null(x)) return()
    sep_cutoffs(vals$re1) <- x
})

# ------------------------------------------------------------------------------
# yield, event and mahal plot
# ------------------------------------------------------------------------------
observe({
    if (is.null(x <- vals$re2)) return()
    output$plot_plotYields <- renderPlot(plotYields(x, input$yp_which))
    output$plot_plotEvents <- renderPlot(plotEvents(x, input$ep_which, as.numeric(input$n_events)))
    output$plot_plotMahal <- renderPlot(plotMahal(x, input$mhl_which, input$input_mhlCofactor))
})

# inputSelect choices for event and mahal plot
observe({
    if (is.null(vals$re2)) return()
    vals$ids <- sort(unique(bc_ids(vals$re2)))
    vals$ep_choices  <- vals$ids
    vals$mhl_choices <- vals$ids[vals$ids != 0]
})

# summary table: IDs | Counts | Cutoffs | Yields
output$table_summary <- DT::renderDataTable({ 
    if (is.null(vals$re2)) return()
    summary_tbl(vals$re2) 
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
