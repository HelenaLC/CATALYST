# ------------------------------------------------------------------------------
# before vs. after compensation
# ------------------------------------------------------------------------------
panel_scatters <- function(samples) {
    fluidPage(
        # sample selection & next / previous sample buttons
        fluidRow(
            align="center",
            div(style=inlineCenter,
                bsButton(
                    inputId="prev_smplComp", 
                    label=NULL,
                    icon=icon("chevron-left"),
                    style="default",
                    size="extra-small")),
            div(style=inline,
                selectInput(
                    inputId="select_smplComp",
                    label=NULL,
                    choices=samples,
                    selected=samples[1],
                    width="320px")),
            div(style=inlineCenter,
                bsButton(
                    inputId="next_smplComp", 
                    label=NULL,
                    icon=icon("chevron-right"),
                    style="default",
                    size="extra-small"))
        ), 
        # channel selection, flip axes button, 
        # cofactor input & view button
        fluidRow(
            align="center",
            div(style=inlineCenter,
                bsButton(
                    inputId="flipAxes", 
                    label=NULL,
                    icon=icon("exchange"),
                    style="primary",
                    size="extra-small")),
            bsTooltip(
                id="flipAxes", 
                title="Swap axes", 
                placement="right", 
                trigger="hover"),
            div(style=inline,
                h5(strong(" X-axis:"), align="right")),
            div(style=inlineTop,
                selectInput(
                    inputId="scatterCh1",
                    label=NULL, 
                    choices=NULL,
                    width="120px")),
            div(style=inline,
                h5(strong(" Y-axis:"), align="right")),
            div(style=inlineTop,
                selectInput(
                    inputId="scatterCh2", 
                    label=NULL, 
                    choices=NULL,
                    width="120px")),
            div(style=inline,
                h5(strong(" Cofactor:"), align="right")),
            div(style=inlineTop,
                numericInput(
                    inputId="cfComp", 
                    label=NULL,
                    value=5, 
                    min=1, 
                    width="91px"))
        ),
        # display current spill & spill adjustment
        fluidRow(
            tags$head(tags$style(type="text/css", "#text_spill {
                height:35px; width:92px; color:red; font-size:14px; padding:6px")), 
            align="center",
            div(style=inline,
                h5(strong("Spillover:"), align="right")),
            div(style=inlineTop,
                verbatimTextOutput("text_spill")),
            div(style=inline,   
                h5(strong( "Enter new:"), align="right")),
            # numericInput: new spill value  
            div(style=inlineTop,
                numericInput(
                    inputId="newSpill", 
                    label=NULL, 
                    value=NULL, 
                    min=0, 
                    max=100, 
                    step=.01, 
                    width="120px")), 
            # bsButton: adjust spill of current interaction
            div(style=inlineCenter,
                bsButton(
                    inputId="adjustSpill", 
                    label=NULL,
                    icon=icon("reply"), 
                    style="success",
                    size="extra-small")),
            bsTooltip(id="adjustSpill", 
                      title="Adjust",
                      placement="right", 
                      trigger="hover"),
            # bsButton: revert current adjustment
            div(style=inlineCenter,
                bsButton(
                    inputId="revert", 
                    label=NULL,
                    icon=icon("reply"),
                    style="warning",
                    size="extra-small")),
            bsTooltip(id="revert", 
                      title="Revert", 
                      placement="right", 
                      trigger="hover"),
            # bsButton: revert all adjustments
            div(style=inlineCenter,
                bsButton(
                    inputId="revertAll", 
                    label=NULL,
                    icon=icon("reply-all"),
                    style="danger",
                    size="extra-small")),
            bsTooltip(id="revertAll", 
                      title="Revert all", 
                      placement="right", 
                      trigger="hover")
        ),
        # before vs. after compensaton scatters
        fluidRow(
            align="center",
            shinydashboard::box(
                width=12, 
                fluidPage( 
                    tags$head(tags$style("#text_info1 {
                        font-size:20px; color:blue; background-color:white")), 
                    tags$head(tags$style("#text_info2 {
                        font-size:20px; color:blue; background-color:white")), 
                    tags$head(tags$style("#compScatter1{
                        float:center; height:100vh !important}")),
                    tags$head(tags$style("#compScatter2{
                        float:center; height:100vh !important}")),
                    fluidRow(
                        align="center",
                        column(
                            width=6, 
                            fluidRow(
                                plotOutput("compScatter1", 
                                brush="rect1",
                                width="100%")),
                            fluidRow(
                                align="left",
                                column(
                                    width=12,
                                    offset=3,
                                    shinydashboard::box(
                                        verbatimTextOutput("text_info1"),
                                        width=6)))),
                        column(
                            width=6, 
                            fluidRow(
                                plotOutput("compScatter2", 
                                brush="rect2", 
                                width="100%")),
                            fluidRow(
                                align="left",
                                column(
                                    width=12,
                                    offset=3,
                                    shinydashboard::box(
                                        verbatimTextOutput("text_info2"),
                                        width=6))))
                    )
                )
            )
        )
    )
}