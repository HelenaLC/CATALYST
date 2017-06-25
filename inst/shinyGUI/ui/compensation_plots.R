inline <- "display:inline-block"
inline36 <- "display:inline-block; line-height:36px"
inlineTop <- "display:inline-block; vertical-align:top"
inlineTop34 <- "display:inline-block; vertical-align:top; height:34px"

# ------------------------------------------------------------------------------
# plotSpillmat()
panel_plotSpillmat <- fluidPage(
        fluidRow(align="center", 
            tags$head(tags$style("#plot_plotSpillmat{height:100vh !important;}")),
            plotOutput("plot_plotSpillmat")))

# ------------------------------------------------------------------------------
# before vs. after compensation
panel_scatters <- function(samples, channels, x, y) {

    fluidPage(
        #tags$head(tags$style(type="text/css", ".selectize-input {vertical-align:middle")),
        #tags$head(tags$style(type="text/css", "#flipAxes   {height:36px;")),
        #tags$head(tags$style(type="text/css", "#cfComp     {height:35px;")),
        #tags$head(tags$style(type="text/css", "#newSpill   {height:35px;")),
        tags$head(tags$style(type="text/css", "#text_spill {height:35px; width:92px; color:firebrick; font-size:14px; padding:6px")), 
        # sample selection & next / previous buttons
        fluidRow(
            align="center",
            div(style="display:inline-block; vertical-align:center",
                bsButton(inputId="prevSmplComp", 
                         label=NULL,
                         icon=icon("chevron-left"),
                         style="default",
                         size="extra-small")),
            div(style="display:inline-block",
                selectInput(
                    inputId="select_sampleComp",
                    label=NULL,
                    choices=samples,
                    selected=samples[1],
                    width="320px")),
            div(style="display:inline-block; vertical-align:middle",
                bsButton(inputId="nextSmplComp", 
                         label=NULL,
                         icon=icon("chevron-right"),
                         style="default",
                         size="extra-small"))
        ), 
        # channel selection, flip axes button, 
        # cofactor input & view button
        fluidRow(
            align="center",
            div(style=inlineTop,
                bsButton(inputId="flipAxes", 
                         label=NULL,
                         icon=icon("exchange"),
                         style="primary")),
            bsTooltip(id="flipAxes", 
                      title="Swap axes", 
                      placement="right", 
                      trigger="hover"),
            div(style="display:inline-block",
                h5(strong(" X-axis:"), align="right")),
            div(style=inlineTop,
                selectInput(inputId="scatterCh1",
                            label=NULL, 
                            choices=channels,
                            selected=channels[x],
                            width="120px")),
            div(style="display:inline-block",
                h5(strong(" Y-axis:"), align="right")),
            div(style=inlineTop,
                selectInput(inputId="scatterCh2", 
                            label=NULL, 
                            choices=channels,
                            selected=channels[y], 
                            width="120px")),
            div(style="display:inline-block",
                h5(strong(" Cofactor:"), align="right")),
            div(style=inlineTop,
                numericInput(inputId="cfComp", 
                             label=NULL,
                             value=5, 
                             min=1, 
                             width="91px"))
        ),
        # display current spill & spill adjustment
        fluidRow(
            align="center",
            div(style="display:inline-block",
                h5(strong("Spillover:"), align="right")),
            div(style=inlineTop,
                verbatimTextOutput("text_spill")),
            div(style="display:inline-block",   
                h5(strong( "Enter new:"), align="right")),
            div(style=inlineTop,
                numericInput(inputId="newSpill", 
                             label=NULL, 
                             value=NULL, 
                             min=0, 
                             max=100, 
                             step=.01, 
                             width="120px")), 
            
            div(style=inlineTop,
                bsButton(inputId="adjustSpill", 
                         label=NULL,
                         icon=icon("reply"), 
                         style="success")),
            bsTooltip(id="adjustSpill", 
                      title="Adjust",
                      placement="right", 
                      trigger="hover"),
            
            div(style=inlineTop,
                bsButton(inputId="revert", 
                         label=NULL,
                         icon=icon("reply"),
                         style="warning")),
            bsTooltip(id="revert", 
                      title="Revert", 
                      placement="right", 
                      trigger="hover"),
            
            div(style=inlineTop,
                bsButton(inputId="revertAll", 
                         label=NULL,
                         icon=icon("reply-all"),
                         style="danger")),
            bsTooltip(id="revertAll", 
                      title="Revert all", 
                      placement="right", 
                      trigger="hover")
        ),
        # before vs. after compensaton scatters
        fluidRow(
            shinydashboard::box(
                width=12, 
                fluidPage( 
                    fluidRow(
                        align="center",
                        tags$head(tags$style("#compScatter1{float:center; height:100vh !important;}")),
                        tags$head(tags$style("#compScatter2{float:center; height:100vh !important;}")),
                        tags$head(tags$style("#text_info1 {width:200px; float:center; text-align:left; padding-left:20px}")),
                        tags$head(tags$style("#text_info2 {width:200px; float:center; text-align:left; padding-left:20px}")),
                        column(width=6, 
                               fluidRow(plotOutput("compScatter1", 
                                                   brush="rect1",
                                                   width="100%")),
                               fluidRow(verbatimTextOutput("text_info1"))),
                        column(width=6, 
                               fluidRow(plotOutput("compScatter2", 
                                                   brush="rect2", 
                                                   width="100%")),
                               fluidRow(verbatimTextOutput("text_info2")))
                    )
                )
            )
        )
    )
}