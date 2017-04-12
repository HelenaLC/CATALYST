# ------------------------------------------------------------------------------
# before vs. after compensation scatters
# ------------------------------------------------------------------------------

# top-left scatter (ff1 uncompensated)
output$plot_scatter1 <- renderPlot({ 
    input$button_view
    plotScatter(
        es=exprs(vals$ff1), 
        x=isolate(input$input_scatterCh1), 
        y=isolate(input$input_scatterCh2),
        cf=isolate(input$input_cofactor),
        n=25e3)
})

output$text_info1 <- renderText({ 
    text_info(
        vals$ff1, 
        isolate(input$input_cofactor), 
        input$rect1, 
        isolate(input$input_scatterCh1), 
        isolate(input$input_scatterCh2))
})   

# top-right scatter (ff1 compensated)
output$plot_scatter2 <- renderPlot({
    if (is.null(vals$cmp1)) return()
    (input$button_view || input$button_newSpill)
    plotScatter(
        es=exprs(isolate(vals$cmp1)),
        x=isolate(input$input_scatterCh1), 
        y=isolate(input$input_scatterCh2),
        cf=isolate(input$input_cofactor),
        n=25e3)
})
output$text_info2 <- renderText({ 
    input$button_view 
    text_info(
        isolate(vals$cmp1), 
        isolate(input$input_cofactor), 
        input$rect2, 
        isolate(input$input_scatterCh1), 
        isolate(input$input_scatterCh2))
})

# bottom-left scatter (ff2 uncompensated)
output$plot_scatter3 <- renderPlot({ 
    input$button_view 
    plotScatter(
        es=exprs(vals$ff2),
        x=isolate(input$input_scatterCh1), 
        y=isolate(input$input_scatterCh2),
        cf=isolate(input$input_cofactor),
        n=25e3)
})
output$text_info3 <- renderText({ 
    text_info(
        vals$ff2, 
        isolate(input$input_cofactor), 
        input$rect3,
        isolate(input$input_scatterCh1), 
        isolate(input$input_scatterCh2))
}) 

# bottom-left scatter (ff2 uncompensated)
output$plot_scatter4 <- renderPlot({ 
    if (is.null(vals$cmp2)) return()
    (input$button_view || input$button_newSpill)
    plotScatter(
        es=exprs(isolate(vals$cmp2)),
        x=isolate(input$input_scatterCh1), 
        y=isolate(input$input_scatterCh2),
        cf=isolate(input$input_cofactor),
        n=25e3)
})
output$text_info4 <- renderText({ 
    input$button_view 
    text_info(
        isolate(vals$cmp2), 
        isolate(input$input_cofactor), 
        input$rect4, 
        isolate(input$input_scatterCh1), 
        isolate(input$input_scatterCh2))
})
