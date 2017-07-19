# ==============================================================================
# debarcoding tab
# ==============================================================================

debarcodingTab <- fluidPage( 
    tags$style(type="text/css", "a{color:steelblue; font-weight:bold}"),
    shinyjs::useShinyjs(),
    sidebarLayout(
        mainPanel(width=9, 
            tabBox(width=12, 
                tabPanel(
                    title=icon("info-circle"), 
                    uiOutput("debarcoding_guide")),
                tabPanel(
                    title="Yield plot", 
                    uiOutput("yieldPlotPanel")),
                tabPanel(
                    title="Event plot", 
                    uiOutput("eventPlotPanel")),
                tabPanel(
                    title="Mahal plot", 
                    uiOutput("mahalPlotPanel")))),
        sidebarPanel(width=3,
            fileInput(
                inputId="fcsDeba", 
                label="Upload FCS", 
                accept=".fcs", 
                width="100%"),
            uiOutput("debarcodingSidebar1"),
            uiOutput("debarcodingSidebar2"),
            uiOutput("debarcodingSidebar3"))
    )
)

# ------------------------------------------------------------------------------
# 1st sidebar
# ------------------------------------------------------------------------------

debarcodingSidebar1 <- tagList(
    checkboxInput(
        inputId="boxUploadCsv", 
        label="Upload barcoding scheme (CSV)", 
        value=TRUE), 
    uiOutput("debaSchemeCsv"),
    checkboxInput(
        inputId="boxSelectBcChs", 
        label="Select single-positive channels"), 
    uiOutput("selectBcChs"),
    bsButton(
        inputId="buttonDebarcode", 
        label=strong("Debarcode"), 
        style="warning",
        size="small",
        block=TRUE)
)

# ------------------------------------------------------------------------------
# 2nd sidebar
# ------------------------------------------------------------------------------

debarcodingSidebar2 <- tagList(
    hr(style="border-color:black"),
    # automated cutoff estimation
    checkboxInput(
        inputId="box_estCutoffs", 
        label="Estimate separation cutoffs", 
        value=TRUE),
    # population-specific cutoffs
    checkboxInput(
        inputId="box_adjustCutoff", 
        label="Adjust population-specific cutoffs"),
    uiOutput("adjustCutoffUI"),
    checkboxInput(
        inputId="box_globalCutoff", 
        label="Enter global separation cutoff"),
    uiOutput("globalCutoffUI"),
    # mahalanobis distance cutoff
    div(style="display:inline-block; vertical-align:middle; width:75%",
        sliderInput(
            inputId="mahalCutoff", 
            label="Mahalanobis distance threshold",
            min=10, 
            max=100, 
            value=30,
            width="100%")),
    div(style=inlineCenter, 
        bsButton(
            inputId="button_mhlCutoff",
            label=NULL,
            icon=icon("share"),
            style="default",
            size="extra-small")),
    bsTooltip(
        id="button_mhlCutoff",
        title="Apply",
        placement="right"),
    hr(style="border-color:black"),
    # checkboxes for output FCS file naming
    checkboxInput(
        inputId="box_IDsAsNms", 
        label="Use sample IDs as files names", 
        value=TRUE),
    checkboxInput(
        inputId="box_upldNms",  
        label="Upload naming sheet (CSV)"),
    fileInput(
        inputId="input_upldNms", 
        label=helpText("This should be a 2 column sheet with
            sample IDs and the desired output filenames."), 
        accept=".csv"),
    # download buttons
    tags$style(type="text/css", "#dwnld_debaFcs {
        display:inline-block; color:white; width:49%; float:left}"),
    tags$style(type="text/css", "#dwnld_debaPlots {
        display:inline-block; color:white; width:49%; float:right}"),
    downloadButton(
        outputId="dwnld_debaFcs", 
        label="FCS files",
        class="btn-success"), 
    downloadButton(
        outputId="dwnld_debaPlots", 
        label="Plots", 
        class="btn-success"))

adjustCutoffUI <- function(dbFrame, choices) {
    tagList(
        div(style="display:inline-block; vertical-align:top; width:25%",
            selectInput(
                inputId="select_adjustCutoff",
                label=NULL,
                choices=choices)),
        div(style="display:inline-block; width:25%",
            numericInput(
                inputId="input_adjustCutoff",
                label=NULL,
                value=sep_cutoffs(dbFrame)[1],
                min=0, max=1, step=.01)),
        div(style=inlineCenter,
            tagList(
                bsButton(
                    inputId="button_adjustCutoff",
                    label=NULL,
                    icon=icon("share"),
                    style="warning",
                    size="extra-small"),
                bsTooltip(
                    id="button_adjustCutoff",
                    title="Adjust",
                    placement="right")))
    )
}

globalCutoffUI <- tagList(
    div(style="display:inline-block; width:25%", 
        numericInput(
            inputId="input_globalCutoff", 
            label=NULL, value=NULL, 
            min=0, max=1, step=.01)),
    div(style=inlineCenter,
        bsButton(
            inputId="button_globalCutoff", 
            label=NULL, 
            icon=icon("share"), 
            style="warning",
            size="extra-small")),
    bsTooltip(
        id="button_globalCutoff",
        title="Apply", 
        placement="right")
)