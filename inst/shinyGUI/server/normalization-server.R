# ------------------------------------------------------------------------------
# bead vs. dna scatters for bead gating
# ------------------------------------------------------------------------------

output$plot_beadScatter1 <- renderPlot({ 
    (input$button_viewGating || input$button_gateBeads)
    if (is.null(vals$beadCols) || is.null(vals$dnaCols)) return()
    plotScatter(
        es=flowCore::exprs(vals$ff_norm),
        x=vals$beadCols[1], 
        y=vals$dnaCols[1],
        cf=isolate(input$input_cfGating), 
        n=isolate(input$input_nGating)) 
})

output$plot_beadScatter2 <- renderPlot({ 
    (input$button_viewGating || input$button_gateBeads)
    if (is.null(vals$beadCols) || is.null(vals$dnaCols)) return()
    plotScatter(
        es=flowCore::exprs(vals$ff_norm),
        x=vals$beadCols[2], 
        y=vals$dnaCols[1],
        cf=isolate(input$input_cfGating), 
        n=isolate(input$input_nGating)) 
})

output$plot_beadScatter3 <- renderPlot({ 
    (input$button_viewGating || input$button_gateBeads)
    if (is.null(vals$beadCols) || is.null(vals$dnaCols)) return()
    plotScatter(
        es=flowCore::exprs(vals$ff_norm),
        x=vals$beadCols[3], 
        y=vals$dnaCols[1],
        cf=isolate(input$input_cfGating), 
        n=isolate(input$input_nGating)) 
})

output$plot_beadScatter4 <- renderPlot({ 
    (input$button_viewGating || input$button_gateBeads)
    if (is.null(vals$beadCols) || is.null(vals$dnaCols)) return()
    plotScatter(
        es=flowCore::exprs(vals$ff_norm),
        x=vals$beadCols[4], 
        y=vals$dnaCols[1],
        cf=isolate(input$input_cfGating), 
        n=isolate(input$input_nGating)) 
})

output$plot_beadScatter5 <- renderPlot({ 
    (input$button_viewGating || input$button_gateBeads)
    if (is.null(vals$beadCols) || is.null(vals$dnaCols)) return()
    plotScatter(
        es=flowCore::exprs(vals$ff_norm),
        x=vals$beadCols[5], 
        y=vals$dnaCols[1],
        cf=isolate(input$input_cfGating), 
        n=isolate(input$input_nGating)) 
})