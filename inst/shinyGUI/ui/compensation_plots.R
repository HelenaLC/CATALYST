# ------------------------------------------------------------------------------
# estTrim()
panel_estTrim <- fluidPage(
        fluidRow(align="center", 
                 plotOutput("plot_estTrim")))

# ------------------------------------------------------------------------------
# plotSpillmat()
panel_plotSpillmat <- fluidPage(
        fluidRow(align="center", 
                 tags$head(tags$style("#plot_plotSpillmat{height:100vh !important;}")),
                 plotOutput("plot_plotSpillmat")))

# ------------------------------------------------------------------------------
# plotScatters()
panel_plotScatter <- fluidPage(
    fluidRow(align="center", 
             plotOutput("plot_plotScatter")))

# ------------------------------------------------------------------------------
# before vs. after compensation
panel_scatters <- function(choices, ch1, ch2) {
    fluidPage(
        tags$head(tags$style(type="text/css", "#text_spill      {height:34px; width:120px; text-align:center")),
        tags$head(tags$style(type="text/css", "#button_view     {height:34px; width:120px; border-width:1.25px; border-color:gold; background-color:lightyellow")),
        tags$head(tags$style(type="text/css", "#button_newSpill {height:34px; width:120px; border-width:1.25px; border-color:gold; background-color:lightyellow")),
        tags$head(tags$style(type="text/css", "#text_info1 {width:200px; float:center; text-align:right; padding-right:25px}")),
        tags$head(tags$style(type="text/css", "#text_info2 {width:200px; float:center; text-align:right; padding-right:25px}")),
        tags$head(tags$style(type="text/css", "#text_info3 {width:200px; float:center; text-align:right; padding-right:25px}")),
        tags$head(tags$style(type="text/css", "#text_info4 {width:200px; float:center; text-align:right; padding-right:25px}")),
        fluidRow(
            column(12, align="center",
                   div(style="display:inline-block",                     h6("X-axis:", align="right")),
                   div(style="display:inline-block; vertical-align:top", selectInput("input_scatterCh1", NULL, choices=choices, selected=ch1, width="120px")),
                   div(style="display:inline-block",                     h6("Y-axis:", align="right")),
                   div(style="display:inline-block; vertical-align:top", selectInput("input_scatterCh2", NULL, choices=choices, selected=ch2, width="120px")),
                   div(style="display:inline-block",                     h6("Cofactor:", align="right")),
                   div(style="display:inline-block",                     numericInput("input_cofactor", NULL, value=10, min=5, width="90px")), 
                   div(style="display:inline-block; vertical-align:top", actionButton("button_view", strong("View")))),
            column(12, align="center",
                   div(style="display:inline-block",                     h6("Spillover:", align="right")),
                   div(style="display:inline-block; vertical-align:top", verbatimTextOutput("text_spill")),
                   div(style="display:inline-block",                     h6("Enter new:", align="right")),
                   div(style="display:inline-block; vertical-align:top", numericInput("input_newSpill", NULL, value=NULL, min=0, max=100, step=.01, width="90px")), 
                   div(style="display:inline-block; vertical-align:top", actionButton("button_newSpill", strong("Adjust"))))),
            
        fluidRow(align="center",
            column(6, fluidRow(plotOutput("plot_scatter1", brush="rect1", width="378px", height="370px")),
                      fluidRow(verbatimTextOutput("text_info1"))), 
            column(6, fluidRow(plotOutput("plot_scatter2", brush="rect2", width="378px", height="370px")),
                      fluidRow(verbatimTextOutput("text_info2"))),
            column(6, fluidRow(plotOutput("plot_scatter3", brush="rect3", width="378px", height="370px")),
                      fluidRow(verbatimTextOutput("text_info3"))),
            column(6, fluidRow(plotOutput("plot_scatter4", brush="rect4", width="378px", height="370px")),
                      fluidRow(verbatimTextOutput("text_info4")))))
}