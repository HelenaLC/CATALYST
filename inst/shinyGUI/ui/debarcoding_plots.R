# ------------------------------------------------------------------------------
# plotYields()
# ------------------------------------------------------------------------------

yieldPlotPanel <- function(x) { 
    fluidPage(
        column(width=8,
            fluidRow(align="center",
                # previous channel button
                div(style=inlineCenter,
                    shinyBS::bsButton(
                        inputId="prev_yieldPlot",
                        label=NULL,
                        icon=icon("chevron-left"),
                        size="extra-small")),
                # which channel selectInput
                div(style=inline,
                    selectInput(
                        inputId="select_yieldPlot", 
                        label=NULL, 
                        width="120px",
                        choices=x)),
                # next channel button
                div(style=inlineCenter,
                    shinyBS::bsButton(
                        inputId="next_yieldPlot",
                        label=NULL,
                        icon=icon("chevron-right"),
                        size="extra-small"))),
            # plotting window
            fluidRow(align="center",
                tags$head(tags$style("#yieldPlot{height:70vh !important;}")),
                plotOutput("yieldPlot"))),
        column(width=4,
            align="center",
            tags$head(tags$style("#table_summary{height:100vh !important;}")),
            DT::dataTableOutput(
                outputId="table_summary")))
}

# ------------------------------------------------------------------------------
# plotEvents()
# ------------------------------------------------------------------------------

eventPlotPanel <- function(x, y) {
    fluidPage(
        fluidRow(align="center",
            # previous channel button
            div(style=inlineCenter,
                shinyBS::bsButton(
                    inputId="prev_eventPlot", 
                    label=NULL,
                    icon=icon("chevron-left"),
                    size="extra-small")),
            # which channel selectInput
            div(style=inline,
                selectInput(
                    inputId="select_eventPlot", 
                    label=NULL,
                    width="120px", 
                    choices=x,
                    selected=y)),
            # next channel button
            div(style=inlineCenter,
                shinyBS::bsButton(
                    inputId="next_eventPlot", 
                    label=NULL,
                    icon=icon("chevron-right"),
                    size="extra-small")),
            # number of events input
            div(style=inline,
                selectInput(
                    inputId="n_events", 
                    label=NULL, 
                    width="120px", 
                    choices=c(1e2,250,1e3,2500,1e4), 
                    selected=1e3))),
        # plotting window
        fluidRow(align="center", 
            tags$head(tags$style("#eventPlot{height:100vh !important;}")),
            plotOutput(outputId="eventPlot")))
}

# ------------------------------------------------------------------------------
# plotMahal()
# ------------------------------------------------------------------------------

mahalPlotPanel <- function(x, y) {
    fluidPage(
        fluidRow(align="center",
            # previous channel button
            div(style=inlineCenter,
                shinyBS::bsButton(
                    inputId="prev_mahalPlot", 
                    label=NULL, 
                    icon=icon("chevron-left"),
                    style="default",
                    size="extra-small")),
            # which channel input
            div(style=inline,
                selectInput(
                    inputId="select_mahalPlot",
                    label=NULL, 
                    choices=x,
                    selected=y,
                    width="120px")),
            # next channel button
            div(style=inlineCenter,
                shinyBS::bsButton(
                    inputId="next_mahalPlot", 
                    label=NULL, 
                    icon=icon("chevron-right"),
                    style="default",
                    size="extra-small")),
            # cofactor numeric input
            div(style=inline, 
                "Cofactor:"),
            div(style=inlineTop,     
                numericInput(
                    inputId="mahalPlotCofactor", 
                    label=NULL, 
                    value=5, 
                    min=1, 
                    width="60px"))), 
        # plotting window
        fluidRow(align="center", 
            tags$head(tags$style("#mahalPlot{height:100vh !important;}")),
            plotOutput(outputId="mahalPlot")))
}
