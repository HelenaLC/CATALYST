# ==============================================================================
# compensation tab
# ==============================================================================

chList <- function(ff, inds, id) {
    chs <- flowCore::colnames(ff)
    n <- sum(inds)
    chList <- lapply(seq_len(n), function(i) {
        textInput(
            inputId=paste0(id, i),
            label=NULL,
            value=chs[inds][i],
            width="100%")
    })
    do.call(tagList, chList)
}

editFCS <- function() {
    modalDialog(
        fluidPage(
            column(
                width=6,
                style="padding:10px",
                uiOutput("chList1")
            ),
            column(
                width=6,
                style="padding:10px",
                uiOutput("chList2")
            )
        ),
        title=strong("FCS file editing"),
        footer=shinyBS::bsButton(inputId="done", label="Done"),
        size="m"
    )
}

compensationTab <- fluidPage(
    shinyjs::extendShinyjs(text=textInputCol),
    tags$style("#plotSpillmat{height:100vh !important; width:100%}"),
    sidebarLayout(
        position="left",
        mainPanel(
            width=9,
            tabBox(
                width=12,
                tabPanel(
                    title=icon("info-circle"), 
                    uiOutput("compensation_guide")),
                tabPanel(
                   title="Yield plot",
                   uiOutput("yieldPlotPanelComp")),
                tabPanel(
                    title="Spillover matrix", 
                    plotOutput("plotSpillmat")),
                tabPanel(
                    title="Before vs. after scatters", 
                    uiOutput("panel_scatters")))),
        sidebarPanel(
            width=3,
            fileInput(
                inputId="fcsComp", 
                label="Upload FCS", 
                multiple=TRUE,
                accept=".fcs"),
            uiOutput("compensationSidebar1"),
            uiOutput("compensationSidebar2"))))

# ------------------------------------------------------------------------------
# 1st sidebar
# ------------------------------------------------------------------------------

compensationSidebar1 <- tagList(
    tags$head(tags$style("#yieldPlotComp{height:70vh !important;}")),
    tags$head(tags$style("#summaryTblComp{height:100vh !important;}")),
    checkboxInput(
        inputId="checkbox_upldSm", 
        label="Upload spillover matrix (CSV)"),
    uiOutput(outputId="uploadSpillMat"),
    checkboxInput(
        inputId="checkbox_estSm", 
        label="Estimate spill from single-stained controls"),
    uiOutput(outputId="uploadControls"),
    uiOutput(outputId="selectSinglePosChs"),
    uiOutput(outputId="debaParsComp"))

# ------------------------------------------------------------------------------
# 2nd sidebar
# ------------------------------------------------------------------------------

compensationSidebar2 <- tagList(
    hr(style="border-color:black"),
    checkboxInput(
        inputId="box_setToZero", 
        label="Should negative values be set to zero?"),
    tags$style(type="text/css", "#dwnld_comped {
        display:inline-block; color:white; width:49%; float:left}"),
    tags$style(type="text/css", "#dwnld_spillMat {
        display:inline-block; color:white; width:49%; float:right}"),
    downloadButton(
        outputId="dwnld_comped", 
        label="Compensated data",
        class="btn-success"), 
    downloadButton(
        outputId="dwnld_spillMat", 
        label="Spillover matrix", 
        class="btn-success"))

# ------------------------------------------------------------------------------

selectSinglePosChsUI <- function(choices) {
    tagList(
        selectInput(
            inputId="singlePosChs", 
            label="Select single-positive channels", 
            choices=choices,
            multiple=TRUE, selectize=FALSE, size=12),
        bsButton(
            inputId="debarcodeComp", 
            label="Debarcode", 
            style="warning",
            size="small",
            block=TRUE)
    )
}

# ------------------------------------------------------------------------------
# before vs. after compensation scatters
# ------------------------------------------------------------------------------
panel_scatters <- function(samples) {
    fluidPage(
        # sample selection & next / previous sample buttons
        fluidRow(
            align="center",
            div(style="display:inline-block; vertical-align:center",
                bsButton(
                    inputId="prev_smplComp", 
                    label=NULL,
                    icon=icon("chevron-left"),
                    style="default",
                    size="extra-small")),
            div(style="display:inline-block",
                selectInput(
                    inputId="select_smplComp",
                    label=NULL,
                    choices=samples,
                    selected=samples[1],
                    width="320px")),
            div(style="display:inline-block; vertical-align:center",
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
            div(style="display:inline-block; vertical-align:center",
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
            div(style="display:inline-block",
                h5(strong(" X-axis:"), align="right")),
            div(style="display:inline-block; vertical-align:top",
                selectInput(
                    inputId="scatterCh1",
                    label=NULL, 
                    choices=NULL,
                    width="120px")),
            div(style="display:inline-block",
                h5(strong(" Y-axis:"), align="right")),
            div(style="display:inline-block; vertical-align:top",
                selectInput(
                    inputId="scatterCh2", 
                    label=NULL, 
                    choices=NULL,
                    width="120px")),
            div(style="display:inline-block",
                h5(strong(" Cofactor:"), align="right")),
            div(style="display:inline-block; vertical-align:top",
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
            div(style="display:inline-block",
                h5(strong("Spillover:"), align="right")),
            div(style="display:inline-block; vertical-align:top",
                verbatimTextOutput("text_spill")),
            div(style="display:inline-block",   
                h5(strong( "Enter new:"), align="right")),
            # numericInput: new spill value  
            div(style="display:inline-block; vertical-align:top",
                numericInput(
                    inputId="newSpill", 
                    label=NULL, 
                    value=NULL, 
                    min=0, 
                    max=100, 
                    step=.01, 
                    width="120px")), 
            # bsButton: adjust spill of current interaction
            div(style="display:inline-block; vertical-align:center",
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
            div(style="display:inline-block; vertical-align:center",
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
            div(style="display:inline-block; vertical-align:center",
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