# ------------------------------------------------------------------------------
# plotYields()
# ------------------------------------------------------------------------------

yp_panel <- function(choices) { 
    fluidPage(
        column(width=9,
               fluidRow(align="center",
                        # previous channel button
                        div(style="display:inline-block; vertical-align:center",
                            shinyBS::bsButton(inputId="yp_prev",
                                              label=NULL,
                                              icon=icon("chevron-left"),
                                              size="extra-small")),
                        # which channel input
                        div(style="display:inline-block",
                            selectInput(inputId="yp_which", 
                                        label=NULL, 
                                        width="120px",
                                        choices=choices,
                                        selected=choices[1])),
                        # next channel button
                        div(style="display:inline-block; vertical-align:center",
                            shinyBS::bsButton(inputId="yp_next",
                                              label=NULL,
                                              icon=icon("chevron-right"),
                                              size="extra-small"))),
               # plotting window
               fluidRow(align="center",
                        tags$head(tags$style("#plot_plotYields{height:70vh !important;}")),
                        plotOutput("plot_plotYields", width="100%"))),
        column(width=3,
               align="center",
               tags$head(tags$style("#table_summary{width:100%; height:100vh !important;}")),
               #tags$style(HTML("#table_summary {background-color:white}")),
               DT::dataTableOutput("table_summary")))
}

# ------------------------------------------------------------------------------
# plotEvents()
# ------------------------------------------------------------------------------

ep_panel <- function(choices) {
    fluidPage(
        fluidRow(align="center",
                 # previous channel button
                 div(style="display:inline-block; vertical-align:center",
                     shinyBS::bsButton(inputId="ep_prev", 
                                       label=NULL,
                                       icon=icon("chevron-left"),
                                       size="extra-small")),
                 # which channel input
                 div(style="display:inline-block",
                     selectInput(inputId="ep_which", 
                                 label=NULL,
                                 width="120px", 
                                 choices=choices, 
                                 selected=choices[1])),
                 # next channel button
                 div(style="display:inline-block; vertical-align:center",
                     shinyBS::bsButton(inputId="ep_next", 
                                       label=NULL,
                                       icon=icon("chevron-right"),
                                       size="extra-small")),
                 # number of events input
                 div(style="display:inline-block",
                     selectInput(inputId="n_events", 
                                 label=NULL, 
                                 width="120px", 
                                 choices=c(1e2,250,1e3,2500,1e4), 
                                 selected=1e3))),
        # plotting window
        fluidRow(align="center", 
                 tags$head(tags$style("#plot_plotEvents{height:100vh !important;}")),
                 plotOutput("plot_plotEvents", width="100%")))
}

# ------------------------------------------------------------------------------
# plotMahal()
# ------------------------------------------------------------------------------

mhl_panel <- function(choices) {
    fluidPage(
        fluidRow(align="center",
                 # previous channel button
                 div(style="display:inline-block; height:34px",
                     actionButton("mhl_prev", NULL, width="30px", icon("chevron-left"))),
                 # which channel input
                 div(style="display:inline-block; vertical-align:top",
                     selectInput("mhl_which", NULL, width="120px", choices=choices, selected=choices[1])),
                 # next channel button
                 div(style="display:inline-block; height:34px",
                     actionButton("mhl_next", NULL, width="30px", icon("chevron-right"))),
                 # cofactor numeric input
                 div(style="display:inline-block", "Cofactor:"),
                 div(style="display:inline-block; vertical-align:top",     
                     numericInput("input_mhlCofactor", NULL, value=50, min=5, width="60px"))), 
        # plotting window
        fluidRow(align="center", 
                 tags$head(tags$style("#plot_plotMahal{height:100vh !important;}")),
                 plotOutput("plot_plotMahal", width="100%")))
}
