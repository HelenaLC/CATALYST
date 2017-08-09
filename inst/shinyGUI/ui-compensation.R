# ==============================================================================
# compensation tab
# ==============================================================================

editMs <- function(id) {
    # id := 'SS' or 'MP'
    modalDialog(
        fluidPage(
                column(
                    offset=2,
                    width=4,
                    align="center",
                    uiOutput(paste0("duplicateMasses", id))),
                column(
                    width=4,
                    align="center",
                    uiOutput(paste0("duplicateMetals", id)))),
        title=HTML(paste(strong("Duplicate masses detected"), 
            "Please select which channels to keep.", sep="<br/>")),
        footer=shinyBS::bsButton(
            inputId=paste0("msChecked", id), 
            label="Done"),
        size="m")
}

editMets <- function(id, header1, header2) {
    # id := 'FLvsSS' or 'SSvsMP'
    # header1 := 'Fulidigm' or 'Single-stains'
    # header2 := 'Single-stains' or 'Multiplexed'
    modalDialog(
        wellPanel(
            style="background-color:white; border:0px; 
                overflow-y:scroll; max-height:100vh",
            fluidPage(
                fluidRow(
                    column(
                        width=4,
                        p(strong("Mass"))),
                    column(
                        width=4,
                        p(strong(header1))),
                    column(
                        width=4,
                        p(strong(header2)))),
                fluidRow(
                    column(
                        width=4,
                        align="center",
                        uiOutput(paste0("masses", id))),
                    column(
                        width=3,
                        align="center",
                        uiOutput(paste0("metals1", id))),
                    column(
                        width=1,
                        uiOutput(paste0("boxes1", id))),
                    column(
                        width=3,
                        align="center",
                        uiOutput(paste0("metals2", id))),
                    column(
                        width=1,
                        uiOutput(paste0("boxes2", id)))))),
        title=HTML(paste(sep="<br/>", strong("Metal mismatches detected"), 
            "The metals selected here will effect spillover estimation.
            Please check carefully and comfirm which metals were used.")),
        footer=shinyBS::bsButton(
            inputId=paste0("metsChecked", id), 
            label="Done"),
        size="m")
}

msList <- function(n, id) {
    # n := number of channels
    # id := 'FLvsSS' or 'SSvsMP'
    l <- lapply(seq_len(n), function(i) {
        div(style="height:40px; width:100%", 
            verbatimTextOutput(
                output=paste0("mass", id, i)))
    })
    do.call(tagList, l)
}
selectMs <- function(mets, inds, id) {
    l <- lapply(seq_along(inds), function(i) {
        selectInput(
            inputId=paste0("metalSelection", id, i),
            label=NULL,
            choices=mets[inds[[i]]],
            width="100%")
    })
    do.call(tagList, l)
}
metsList <- function(n, id1, id2) {
    # n := number of channels
    # id1 := 'FLvsSS' or 'SSvsMP'
    # id2 := 'ref' or 'dat'
    l <- lapply(seq_len(n), function(i) {
        div(style="height:40px; width:100%", 
            verbatimTextOutput(
                output=paste0("metal", id1, id2, i)))
    })
    do.call(tagList, l)
}
boxList <- function(n, id1, id2) {
    # n := number of channels
    # id1 := 'FLvsSS' or 'SSvsMP'
    # id2 := 'ref' or 'dat'
    val <- id2 == "dat"
    l <- lapply(seq_len(n), function(i) {
        div(style="height:30px; margin:0; padding:0",
            checkboxInput(
                inputId=paste0("box", id1, id2, i),
                label=NULL,
                value=val))
    })
    do.call(tagList, l)
}

compensationTab <- fluidPage(
    shinyjs::extendShinyjs(text=restyleMetals),
    tags$style("#plotSpillmat{height:100vh !important; width:100%}"),
    tags$style("#yieldPlotComp{height:70vh !important}"),
    tags$style("#summaryTblComp{height:100vh !important}"),
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
            checkboxInput(
                inputId="checkbox_upldSm", 
                label="Upload spillover matrix (CSV)"),
            uiOutput(outputId="uploadSpillMat"),
            checkboxInput(
                inputId="checkbox_estSm", 
                label="Estimate spill from controls"),
            uiOutput(outputId="uploadSs"),
            uiOutput(outputId="uploadMp"),
            uiOutput(outputId="selectSinglePosChs"),
            uiOutput(outputId="debaParsComp"),
            uiOutput("compensationSidebar"))))

# ------------------------------------------------------------------------------
# sidebar extension
# ------------------------------------------------------------------------------

compensationSidebar <- tagList(
    hr(style="border-color:black"),
    checkboxInput(
        inputId="box_setToZero", 
        label="Should negative values be set to zero?"),
    tags$style(type="text/css", "#dwnld_comped {
        display:inline-block; color:white; width:35%; float:left}"),
    tags$style(type="text/css", "#dwnld_spillMat {
        display:inline-block; color:white; width:35%; margin-left:5px}"),
    downloadButton(
        outputId="dwnld_comped", 
        label="Compensated data",
        class="btn-success"), 
    downloadButton(
        outputId="dwnld_spillMat", 
        label="Spillover matrix", 
        class="btn-success"),
    # bsButton: "Go to debarcoding"
    div(style=" display:inline-block; width:25%; float:right",
        bsButton(
            inputId="goToDeba",
            label="Go to debarcoding",
            width="100%"))

)

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