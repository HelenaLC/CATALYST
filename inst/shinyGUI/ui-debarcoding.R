# ==============================================================================
# debarcoding tab
# ==============================================================================

debarcodingTab <- fluidPage( 
    tags$style(type="text/css", "a{color:steelblue; font-weight:bold}"),
    sidebarLayout(
        mainPanel(
            width=9, 
            tabBox(
                width=12, 
                tabPanel(
                    title="Yield plot", 
                    uiOutput("yieldPlotPanelDeba")),
                tabPanel(
                    title="Event plot", 
                    uiOutput("eventPlotPanel")),
                tabPanel(
                    title="Mahal plot", 
                    uiOutput("mahalPlotPanel")))),
        sidebarPanel(
            width=3,
            fileInput(
                inputId="fcsDeba", 
                label="Upload FCS", 
                accept=".fcs", 
                width="100%"),
            uiOutput("debarcodingSidebar1"),
            uiOutput("debarcodingSidebar2"))
    )
)

# ------------------------------------------------------------------------------
# 1st sidebar extension
# ------------------------------------------------------------------------------

debarcodingSidebar1 <- tagList(
    fileInput(
        inputId="debaSchemeCsv", 
        label="Upload CSV", 
        accept=".csv"),
    shinyBS::bsButton(
        inputId="debarcodeDeba", 
        label="Debarcode", 
        style="warning",
        size="small",
        block=TRUE,
        disabled=TRUE))

# ------------------------------------------------------------------------------
# 2nd sidebar extension
# ------------------------------------------------------------------------------

debarcodingSidebar2 <- tagList(
    tags$head(tags$style("#yieldPlotDeba{height:70vh !important;}")),
    tags$head(tags$style("#summaryTblDeba{height:100vh !important;}")),
    tags$style(type="text/css", "#dwnld_debaFcs {
        display:inline-block; color:white; width:49%; float:left}"),
    tags$style(type="text/css", "#dwnld_debaPlots {
        display:inline-block; color:white; width:49%; float:right}"),
    debaParsModule(module="Deba"),
    hr(style="border-color:black"),
    # checkboxes for output FCS file naming
    checkboxInput(
        inputId="box_IDsAsNms", 
        label="Use sample IDs as files names", 
        value=TRUE),
    checkboxInput(
        inputId="box_upldNms",  
        label="Upload naming sheet (CSV)"),
    # fileInput for naming sheet CSV
    fileInput(
        inputId="input_upldNms", 
        label=helpText("This should be a 2 column sheet with
            sample IDs and the desired output filenames."), 
        accept=".csv"),
    # downloadButton
    downloadButton(
        outputId="dwnld_debaFcs", 
        label="FCS files", 
        class="btn-success",
        width="100%"),
    # downloadButton
    downloadButton(
        outputId="dwnld_debaPlots", 
        label="Plots", 
        class="btn-success",
        width="100%"))

# ------------------------------------------------------------------------------
# plotEvents()
# ------------------------------------------------------------------------------

eventPlotPanel <- function(x, y) {
    fluidPage(
        fluidRow(
            align="center",
            # previous channel button
            div(style="display:inline-block; vertical-align:center",
                shinyBS::bsButton(
                    inputId="prev_eventPlot", 
                    label=NULL,
                    icon=icon("chevron-left"),
                    size="extra-small")),
            # which channel selectInput
            div(style="display:inline-block",
                selectInput(
                    inputId="select_eventPlot", 
                    label=NULL,
                    width="120px", 
                    choices=x,
                    selected=y)),
            # next channel button
            div(style="display:inline-block; vertical-align:center",
                shinyBS::bsButton(
                    inputId="next_eventPlot", 
                    label=NULL,
                    icon=icon("chevron-right"),
                    size="extra-small")),
            # number of events input
            div(style="display:inline-block",
                selectInput(
                    inputId="n_events", 
                    label=NULL, 
                    width="120px", 
                    choices=c(1e2,250,1e3,2500,1e4), 
                    selected=1e3))),
        # plotting window
        fluidRow(
            align="center", 
            tags$head(tags$style("#eventPlot{height:100vh !important;}")),
            plotOutput(outputId="eventPlot")))
}

# ------------------------------------------------------------------------------
# plotMahal()
# ------------------------------------------------------------------------------

mahalPlotPanel <- function(x, y) {
    fluidPage(
        fluidRow(
            align="center",
            # previous channel button
            div(style="display:inline-block; vertical-align:center",
                shinyBS::bsButton(
                    inputId="prev_mahalPlot", 
                    label=NULL, 
                    icon=icon("chevron-left"),
                    style="default",
                    size="extra-small")),
            # which channel input
            div(style="display:inline-block",
                selectInput(
                    inputId="select_mahalPlot",
                    label=NULL, 
                    choices=x,
                    selected=y,
                    width="120px")),
            # next channel button
            div(style="display:inline-block; vertical-align:center",
                shinyBS::bsButton(
                    inputId="next_mahalPlot", 
                    label=NULL, 
                    icon=icon("chevron-right"),
                    style="default",
                    size="extra-small")),
            # cofactor numeric input
            div(style="display:inline-block", 
                "Cofactor:"),
            div(style="display:inline-block; vertical-align:top",     
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
