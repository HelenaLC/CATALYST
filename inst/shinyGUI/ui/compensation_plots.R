inline <- "display:inline-block"
inlineTop <- "display:inline-block; vertical-align:top"

# ------------------------------------------------------------------------------
# estTrim()
panel_estTrim <- fluidPage(
    tags$head(tags$style(type="text/css", "#button_estTrim {
        height:34px; width:90px; border-width:1.25px;
        border-color:gold; background-color:lightyellow")),
    fluidRow(align="center",
        # input estTrim min
        div(style=inline, h6("From:")),
        div(style=inline, 
            numericInput("estTrim_min", NULL, width="90px",
                value=.04, min=0, max=.5, step=.01)),
        # input estTrim max
        div(style=inline, h6("To:")),
        div(style=inline, 
            numericInput("estTrim_max", NULL, width="90px",
                value=.18, min=.01, max=.5, step=.01)),
        # input estTrim step
        div(style=inline, h6("Step:")),
        div(style=inline, 
            numericInput("estTrim_step", NULL, width="90px", 
                value=.02, min=.01, max=.5, step=.01)),
        # go button
        div(style=inlineTop, 
            actionButton("button_estTrim", strong("Go")))),
    fluidRow(align="center",
        tags$head(tags$style("#plot_estTrim{height:100vh !important;}")),
        # plotting window
        plotlyOutput("plot_estTrim", width="100%")
    )
)

# ------------------------------------------------------------------------------
# plotSpillmat()
panel_plotSpillmat <- fluidPage(
        fluidRow(align="center", 
            tags$head(tags$style("#plot_plotSpillmat{height:100vh !important;}")),
            plotOutput("plot_plotSpillmat")))

# ------------------------------------------------------------------------------
# before vs. after compensation
panel_scatters <- function(choices, ch1, ch2) {
    fluidPage(
        tags$head(tags$style(type="text/css", "#flipAxes        {height:34px;")),
        tags$head(tags$style(type="text/css", "#button_view     {height:34px;")),
        tags$head(tags$style(type="text/css", "#text_spill      {height:34px; width:120px; text-align:center")),
        tags$head(tags$style(type="text/css", "#button_newSpill {height:34px; width:120px; border-width:1.25px; border-color:gold; background-color:lightyellow")),
        tags$head(tags$style(type="text/css", "#text_info1 {width:200px; float:center; text-align:right; padding-right:25px}")),
        tags$head(tags$style(type="text/css", "#text_info2 {width:200px; float:center; text-align:right; padding-right:25px}")),
        tags$head(tags$style(type="text/css", "#text_info3 {width:200px; float:center; text-align:right; padding-right:25px}")),
        tags$head(tags$style(type="text/css", "#text_info4 {width:200px; float:center; text-align:right; padding-right:25px}")),
        fluidRow(
            column(12, align="center",
                   div(style="display:inline-block; vertical-align:top", 
                       actionButton("flipAxes", NULL, icon("exchange"))),
                   div(style=inline,    h5(strong("X-axis:"), align="right")),
                   div(style=inlineTop, selectInput("input_scatterCh1", NULL, choices=choices, selected=ch1, width="120px")),
                   div(style=inline,    h5(strong("Y-axis:"), align="right")),
                   div(style=inlineTop, selectInput("input_scatterCh2", NULL, choices=choices, selected=ch2, width="120px")),
                   div(style=inline,    h5(strong("Cofactor:"), align="right")),
                   div(style=inline,    numericInput("input_cofactor", NULL, value=10, min=5, width="60px")), 
                   div(style="display:inline-block; vertical-align:top", 
                       actionButton("button_view", NULL, icon("eye")))),
            column(12, align="center",
                   div(style=inlineTop, h6("Spillover:", align="right")),
                   div(style=inlineTop, verbatimTextOutput("text_spill")),
                   div(style=inline,    h6("Enter new:", align="right")),
                   div(style=inlineTop, numericInput("input_newSpill", NULL, value=NULL, min=0, max=100, step=.01, width="90px")), 
                   div(style=inlineTop, actionButton("button_newSpill", strong("Adjust"))))),
            
        fluidRow(align="center",
            tags$head(tags$style("#plot_scatter1{float:right;}")),
            tags$head(tags$style("#plot_scatter3{float:right;}")),
            tags$head(tags$style("#plot_scatter2{float:left}")),
            tags$head(tags$style("#plot_scatter4{float:left}")),
            column(6, fluidRow(plotOutput("plot_scatter1", brush="rect1", width="400px", height="400px")),
                      fluidRow(verbatimTextOutput("text_info1"))), 
            column(6, fluidRow(plotOutput("plot_scatter2", brush="rect2", width="400px", height="400px")),
                      fluidRow(verbatimTextOutput("text_info2"))),
            column(6, fluidRow(plotOutput("plot_scatter3", brush="rect3", width="405px", height="400px")),
                      fluidRow(verbatimTextOutput("text_info3"))),
            column(6, fluidRow(plotOutput("plot_scatter4", brush="rect4", width="378px", height="370px")),
                      fluidRow(verbatimTextOutput("text_info4")))))
}