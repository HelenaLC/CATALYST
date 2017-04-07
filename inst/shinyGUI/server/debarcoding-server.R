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

output$select_bcChs <- renderUI({
    x <- input$box_bcChs
    if (is.null(x) || x == 0) return()
        selectInput("input_bcChs", NULL, choices=flowCore::colnames(vals$ff1), 
            multiple=TRUE, selectize=FALSE, size=12)
})

observe({
    x <- input$box_estCutoffs
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_adjustCutoff", value=FALSE)
    updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
    vals$re2 <- CATALYST::estCutoffs(x=vals$re1)
})

observe({
    x <- input$box_adjustCutoff
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_estCutoffs",   value=FALSE)
    updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
    # synchronize selectInput w/ yield plot
    x <- input$yp_which
    if (is.null(x) || x == 0) return()
    updateSelectInput(session, "select_adjustCutoff", selected=input$yp_which)
})

observe({
    x <- input$box_globalCutoff
    if (is.null(x) || x == 0) return()
    updateCheckboxInput(session, "box_estCutoffs",   value=FALSE)
    updateCheckboxInput(session, "box_adjustCutoff", value=FALSE)
})

# ------------------------------------------------------------------------------
# event and mahal plot
# ------------------------------------------------------------------------------

output$plot_plotEvents <- renderPlot({
    if (is.null(vals$re3) 
        || is.null(input$ep_which) 
        || is.null(input$n_events))
        return()
    CATALYST::plotEvents(
        x=vals$re3, 
        which=input$ep_which,
        n_events=as.numeric(input$n_events))
})

output$plot_plotMahal <- renderPlot({
    if (is.null(vals$re3) 
        || is.null(input$mhl_which) 
        || is.null(input$input_mhlCofactor))
        return()
    CATALYST::plotMahal(
        x=vals$re3, 
        which=input$mhl_which,
        cofactor=input$input_mhlCofactor)
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
    if (is.null(input$yp_which) 
        || input$yp_which == 0) return() 
    updateSelectInput(session, "ep_which",  selected=input$yp_which)
    updateSelectInput(session, "mhl_which", selected=input$yp_which)
})

observe({
    if (is.null(input$ep_which) || input$ep_which == 0) return() 
    updateSelectInput(session, "yp_which",  selected=input$ep_which)
    updateSelectInput(session, "mhl_which", selected=input$ep_which)
})

observe({
    if (is.null(input$mhl_which)) return()
    updateSelectInput(session, "yp_which",  selected=input$mhl_which)
    updateSelectInput(session, "ep_which",  selected=input$mhl_which)
})

# inputSelect choices
observe({
    if (is.null(vals$re3)) return()
    vals$ids <- sort(unique(bc_ids(vals$re3)))
    vals$ep_choices  <- vals$ids
    vals$mhl_choices <- vals$ids[vals$ids != 0]
})

# summary table: IDs | Counts | Cutoffs | Yields
output$table_summary <- DT::renderDataTable({ 
    if (is.null(vals$re3)) return()
    summary_tbl(vals$re3) 
})
