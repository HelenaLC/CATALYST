# ------------------------------------------------------------------------------
# plotYields()
# ------------------------------------------------------------------------------

yp_panel <- function(choices) { 
    fluidPage(
        column(width=9,
            fluidRow(align="center",
                # previous channel button
                div(style=inlineCenter,
                    shinyBS::bsButton(
                        inputId="yp_prev",
                        label=NULL,
                        icon=icon("chevron-left"),
                        size="extra-small")),
                # which channel selectInput
                div(style=inline,
                    selectInput(
                        inputId="yp_which", 
                        label=NULL, 
                        width="120px",
                        choices=choices,
                        selected=choices[1])),
                # next channel button
                div(style=inlineCenter,
                    shinyBS::bsButton(
                        inputId="yp_next",
                        label=NULL,
                        icon=icon("chevron-right"),
                        size="extra-small"))),
            # plotting window
            fluidRow(align="center",
                tags$head(tags$style("#plotYields{height:70vh !important;}")),
                plotOutput("plotYields"))),
        column(width=3,
            align="center",
            tags$head(tags$style("#table_summary{height:100vh !important;}")),
            DT::dataTableOutput(
                outputId="table_summary")))
}

# ------------------------------------------------------------------------------
# plotEvents()
# ------------------------------------------------------------------------------

ep_panel <- function(choices) {
    fluidPage(
        fluidRow(align="center",
            # previous channel button
            div(style=inlineCenter,
                shinyBS::bsButton(
                    inputId="ep_prev", 
                    label=NULL,
                    icon=icon("chevron-left"),
                    size="extra-small")),
            # which channel selectInput
            div(style=inline,
                selectInput(
                    inputId="ep_which", 
                    label=NULL,
                    width="120px", 
                    choices=choices, 
                    selected=choices[1])),
            # next channel button
            div(style=inlineCenter,
                shinyBS::bsButton(
                    inputId="ep_next", 
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
            tags$head(tags$style("#plotEvents{height:100vh !important;}")),
            plotOutput(outputId="plotEvents")))
}

# ------------------------------------------------------------------------------
# plotMahal()
# ------------------------------------------------------------------------------

mhl_panel <- function(choices) {
    fluidPage(
        fluidRow(align="center",
            # previous channel button
            div(style=inlineCenter,
                shinyBS::bsButton(
                    inputId="mhl_prev", 
                    label=NULL, 
                    icon=icon("chevron-left"),
                    style="default",
                    size="extra-small")),
            # which channel input
            div(style=inline,
                selectInput(
                    inputId="mhl_which",
                    label=NULL, 
                    choices=choices, 
                    selected=choices[1], 
                    width="120px")),
            # next channel button
            div(style=inlineCenter,
                shinyBS::bsButton(
                    inputId="mhl_next", 
                    label=NULL, 
                    icon=icon("chevron-right"),
                    style="default",
                    size="extra-small")),
            # cofactor numeric input
            div(style=inline, 
                "Cofactor:"),
            div(style=inlineTop,     
                numericInput(
                    inputId="input_mhlCofactor", 
                    label=NULL, 
                    value=50, 
                    min=5, 
                    width="60px"))), 
        # plotting window
        fluidRow(align="center", 
            tags$head(tags$style("#plotMahal{height:100vh !important;}")),
            plotOutput(outputId="plotMahal")))
}
