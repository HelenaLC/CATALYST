inline <- "display:inline-block"
inlineTop <- "display:inline-block; vertical-align:top"

# ------------------------------------------------------------------------------
# bead gating
# ------------------------------------------------------------------------------
gating_panel <- fluidPage(
    tags$head(tags$style(type="text/css", "#button_viewGating {height:34px; width:120px; border-width:1.25px; border-color:gold; background-color:lightyellow")),
    fluidRow(
        column(12, align="center",
            div(style=inline, h6(strong("Cofactor:"), align="right")),
            div(style=inline, numericInput("input_cfGating", NULL, value=5, min=1, width="90px")),
            div(style=inline, h6(strong("# Events:"), align="right")),
            div(style=inline, numericInput("input_nGating", NULL, value=25e3, min=1e4, width="90px")),
            div(style=inlineTop, actionButton("button_viewGating", strong("View"))))),
    fluidRow(align="center",
        column(2, plotOutput("plot_beadScatter1", brush="gate1", width="200px", height="200px")),
        column(2, plotOutput("plot_beadScatter2", brush="gate2", width="200px", height="200px")),
        column(2, plotOutput("plot_beadScatter3", brush="gate3", width="200px", height="200px")),
        column(2, plotOutput("plot_beadScatter4", brush="gate4", width="200px", height="200px")),
        column(2, plotOutput("plot_beadScatter5", brush="gate5", width="200px", height="200px"))))


