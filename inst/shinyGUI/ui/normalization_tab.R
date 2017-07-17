# ==============================================================================
# normalization tab
# ==============================================================================

normalizationTab <- fluidPage(
    useShinyjs(),
    extendShinyjs(text=collapseBox),
    tags$head(tags$style("#plot_smoothedBeads {height:100vh !important;}")),
    tags$head(tags$style(HTML(".small-box{height:96px; margin-bottom:0px}"))),
    tags$head(tags$style("#dwnld_normResults {color:white; width:100%}")),
    tags$style(HTML(".shiny-plot-output {display:inline-block}")),
    fluidRow(
        column(width=3, 
            style="padding:0px", 
            shinydashboard::box(
                title="Upload FCS",
                solidHeader=TRUE,
                status="warning",
                id="box_1",
                width=12, 
                collapsible=TRUE,
                fileInput(
                    inputId="fcsNorm", 
                    label=NULL,
                    multiple=TRUE,
                    accept=".fcs"))),
        column(width=3, 
            style="padding:0px", 
            uiOutput("box2")),
        column(width=3, 
            style="padding:0px", 
            uiOutput("box3")),
        column(width=3, 
            style="padding:0px", 
            uiOutput("box4"))),
    fluidRow(
        uiOutput("box_beadGating")),
    fluidRow(
        uiOutput("box_smoothedBeads")),
    fluidRow(
        uiOutput("box_beadRemoval"))
)

# box 2: select data to normalize to
box2 <- shinydashboard::box(
    title="Normalization baseline",
    solidHeader=TRUE,
    status="warning",
    id="box_2",
    width=12,
    collapsible=TRUE,
    checkboxInput(
        inputId="box_normToCurrent", 
        label="Normalize to median level of current files"),
    checkboxInput(
        inputId="box_uploadNormTo", 
        label="Upload FCS file(s) of beads to normalize to"),
    uiOutput("input_NormTo")
)

# box 3: bead selection
box3 <- shinydashboard::box(
    title="Bead selection",
    solidHeader=TRUE,
    status="warning",
    id="box_3",
    width=12,
    collapsible=TRUE,
    selectInput(
        inputId="select_beads", 
        label=NULL,
        choices=c("",
            "DVS Beads (140, 151, 153, 165, 175)"="dvs",
            "Beta Beads (139, 141, 159, 169, 175)"="beta",
            "Custom"="custom")),
    uiOutput("select_customBeads")
)

# selectInput for custom beads
selectInput_customBeads <- function(ff) {
    pars <- flowCore::colnames(ff)
    desc <- flowCore::parameters(ff)$desc
    fluidRow(
        column(
            width=12,
            div(style="display:inline-block; width:80%", 
                selectizeInput(
                    inputId="input_customBeads", 
                    label=NULL, 
                    multiple=TRUE, 
                    choices=setNames(
                        object=as.list(pars), 
                        nm=paste0(desc, " [", pars, "]")))),
            div(style="display:inline-block; vertical-align:top; width:19%", 
                actionButton(
                    inputId="button_customBeads", 
                    label=NULL, 
                    icon=icon("share")))
        )
    )
}

# box 4: gating yield and download results
box4 <- shinydashboard::box(
    id="box_4",
    width=12, 
    collapsible=TRUE,
    div(style="display:inline-block; float:right; width:49%; margin-bottom:5px",
        valueBox(
            value=textOutput("howManyGated"), 
            subtitle="samples gated",
            icon=icon("hashtag"),
            color="teal",
            width=NULL)),
    div(style="display:inline-block; float:left; width:49%; margin-bottom:5px",
        valueBox(
            value=textOutput("gatingYield"),
            subtitle="gating yield",
            icon=icon("percent"),
            color="teal", 
            width=NULL)),
    checkboxInput(
        inputId="box_removeBeads",
        label="Should beads be removed?"),
    downloadButton(
        outputId="dwnld_normResults", 
        label="Download results", 
        class="btn-success")
)

# ------------------------------------------------------------------------------
# bead gating panel  
# ------------------------------------------------------------------------------
box_beadGating <- function(samples, selected) {
    shinydashboard::box(
        id="beadGating",
        width=12, 
        collapsible=TRUE,
        fluidPage(
            fluidRow(
                align="center",
                # previous sample button
                div(style=inlineCenter,
                    bsButton(
                        inputId="prevSmplGating", 
                        label=NULL,
                        icon=icon("chevron-left"),
                        style="default",
                        size="extra-small")),
                # sample selection
                div(style=inline,  
                    selectInput(
                        inputId="selectSmplGating",
                        label=NULL,
                        choices=samples,
                        selected=samples[[selected]],
                        width="320px")),
                # next sample button
                div(style=inlineCenter,
                    bsButton(
                        inputId="nextSmplGating", 
                        label=NULL,
                        icon=icon("chevron-right"),
                        style="default",
                        size="extra-small")),
                div(style=inlineTop,
                    bsButton(
                        inputId="gateBeads", 
                        label=NULL,
                        icon=icon("object-ungroup"),
                        style="default",
                        disabled=TRUE)),
                bsTooltip(
                    id="gateBeads", 
                    title="Gate",
                    placement="right", 
                    trigger="hover")
            ),
            # beads vs. dna scatters
            fluidRow(
                align="center",
                uiOutput(outputId="beadScatters")
            )
        )
    )
}

# ------------------------------------------------------------------------------
# smoothed beads panel 
# ------------------------------------------------------------------------------
box_smoothedBeads <- 
    shinydashboard::box(
        id="smoothedBeads",
        width=12, 
        collapsible=TRUE, 
        fluidPage(
            tags$style("#beadsVsBeads{height:100vh !important}"),
            plotOutput(outputId="plot_smoothedBeads", width="100%")
        )
    )

mhlCutoffNormUI <- function(samples) {
    tagList(
        # previous sample button
        div(style=inlineCenter,
            bsButton(
                inputId="prevSmplMhl", 
                label=NULL,
                icon=icon("chevron-left"),
                style="default",
                size="extra-small")),
        # sample selection
        div(style=inline,  
            selectInput(
                inputId="selectSmplMhl",
                label=NULL,
                choices=samples,
                width="320px")),
        # next sample button
        div(style=inlineCenter,
            bsButton(
                inputId="nextSmplMhl", 
                label=NULL,
                icon=icon("chevron-right"),
                style="default",
                size="extra-small")),
        div(style="display:inline-block; width:25%; vertical-align:middle",
            sliderInput(
                inputId="mhlCutoffNorm", 
                label=NULL, 
                min=0, 
                max=100, 
                value=0,
                step=1,
                width="100%")),
        div(style=inlineCenter,
            bsButton(
                inputId="applyMhlCutoffNorm", 
                label=NULL,
                icon=icon("share"),
                style="default",
                size="extra-small")),
        bsTooltip(
            id="applyMhlCutoffNorm", 
            title="Apply",
            placement="right", 
            trigger="hover")
)}  

box_beadRemoval <- 
    shinydashboard::box(
        width=12,
        collapsible=TRUE,
        fluidPage(
            fluidRow(
                align="center",
                uiOutput(outputId="mhlCutoffNormUI")
            ),
            fluidRow(
                tags$style("#beadsVsBeads{height:100vh !important}"),
                plotOutput(outputId="beadsVsBeads", width="100%")
            )
        )
    )