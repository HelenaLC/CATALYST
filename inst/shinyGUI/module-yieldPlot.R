yieldPlotModule <- function(sampleIDs, module) {
    # sampleIDs = choices for selectInput
    # module = "Deba" or "Comp"
    fluidPage(
        column(
            width=8,
            fluidRow(
                align="center",
                # previous channel button
                div(style="display:inline-block; vertical-align:middle",
                    bsButton(
                        inputId=paste0("prev_yieldPlot", module),
                        label=NULL,
                        icon=icon("chevron-left"),
                        size="extra-small")),
                # which channel selectInput
                div(style="display:inline-block; width:25%",
                    selectInput(
                        inputId=paste0("select_yieldPlot", module),
                        label=NULL, 
                        width="100%",
                        choices=sampleIDs)),
                # next channel button
                div(style="display:inline-block; vertical-align:middle",
                    bsButton(
                        inputId=paste0("next_yieldPlot", module),
                        label=NULL,
                        icon=icon("chevron-right"),
                        size="extra-small"))),
            # plotting window
            fluidRow(
                align="center",
                plotOutput(outputId=paste0("yieldPlot", module)))
        ),
        # summary table: IDs | Counts | Cutoffs | Yields
        column(
            width=4,
            align="center",
            DT::dataTableOutput(outputId=paste0("summaryTbl", module))
        )
    )
}