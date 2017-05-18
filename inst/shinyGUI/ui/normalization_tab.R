inline      <- "display:inline-block;"
inlineTop   <- "display:inline-block; vertical-align:top"
inlineTop34 <- "display:inline-block; vertical-align:top; height:34px"

collapseBox <- "shinyjs.collapse = function(boxId) {
$('#'+boxId).closest('.box').find('[data-widget=collapse]').click();}"

# ==============================================================================
# normalization tab
# ==============================================================================

normalizationTab <- fluidPage(
    shinyjs::useShinyjs(),
    shinyjs::extendShinyjs(text=collapseBox),
    # ----- styles -----
    tags$head(tags$style(HTML(".small-box{height:96px; 
                              margin-bottom:0px}"))),
    tags$head(tags$style("#button_customBeads{height:34px}")),
    tags$head(tags$style("#plot_beadScatter1{float:left}")),
    tags$head(tags$style("#plot_beadScatter2{float:left}")),
    tags$head(tags$style("#plot_beadScatter3{float:left}")),
    tags$head(tags$style("#plot_beadScatter4{float:left}")),
    tags$head(tags$style("#plot_beadScatter5{float:left}")),
    tags$head(tags$style("#plot_smoothedBeads{height:100vh !important;}")),
    tags$head(tags$style("#dwnld_normResults {color:white; width:100%}")),
    # -----
    fluidRow(
        column(3, 
               style="padding:0px", 
               shinydashboard::box(
                   id="box_1",
                   width=12, 
                   collapsible=TRUE,
                   fileInput(
                       inputId="fcsNorm", 
                       label="Upload FCS", 
                       multiple=TRUE,
                       accept=".fcs"))),
        column(3, 
               style="padding:0px", 
               uiOutput("box2")),
        column(3, 
               style="padding:0px", 
               uiOutput("box3")),
        column(3, 
               style="padding:0px", 
               uiOutput("box4"))
    ),
    fluidRow(
        uiOutput("box_beadGating")),
    fluidRow(
        uiOutput("box_smoothedBeads"))
    )

# box 2: select data to normalize to -------------------------------------------
box2 <- shinydashboard::box(
    id="box_2",
    width=12,
    collapsible=TRUE,
    checkboxInput(
        inputId="box_NormToCurrent", 
        label="Normalize to median level of current files"),
    checkboxInput(
        inputId="box_uploadNormTo", 
        label="Upload FCS file(s) of beads to normalize to"),
    uiOutput("input_NormTo")
)

# box 3: bead selection --------------------------------------------------------
box3 <- shinydashboard::box(
    id="box_3",
    width=12,
    collapsible=TRUE,
    selectInput(
        inputId="select_beads", 
        label=h5(strong("Select beads")), 
        choices=c("", 
                  "DVS Beads (140, 151, 153, 165, 175)" = "dvs",
                  "Beta Beads (139, 141, 159, 169, 175)" = "beta",
                  "Custom" = "custom")),
    uiOutput("select_customBeads")
)

# selectInput for custom beads
selectInput_customBeads <- function(x) {
    nms <- flowCore::colnames(x)
    des <- flowCore::parameters(x)$desc
    fluidRow(
        column(12,
               div(style="display:inline-block; width:80%", 
                   selectizeInput(
                       inputId="input_customBeads", 
                       label=NULL, 
                       multiple=TRUE, 
                       choices=setNames(
                           object=as.list(nms), 
                           nm=paste0(des, " [", nms, "]")))),
               div(style="display:inline-block; vertical-align:top; width:19%", 
                   actionButton(
                       inputId="button_customBeads", 
                       label="Gate", 
                       icon=icon("object-ungroup"), 
                       width="100%"))
        ))
}

# box 4: gating yield and download results -------------------------------------
box4 <- shinydashboard::box(
    id="box_4",
    width=12, 
    collapsible=TRUE,
    div(style="display:inline-block; float:right; width:49%; margin-bottom:5px",
        valueBox(
            value=textOutput("howManyGated"), 
            subtitle="samples gated",
            icon=icon("hashtag"),
            color="light-blue",
            width=NULL)),
    div(style="display:inline-block; float:left; width:49%; margin-bottom:5px",
        valueBox(
            value=textOutput("gatingYield"),
            subtitle="gating yield",
            icon=icon("percent"),
            color="light-blue", 
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
        width=12, 
        collapsible=TRUE,
        footer=em(strong("Brush to gate beads."), 
                  "Ideally, bead-bead and cell-bead doublets should be excluded.
                  However, if bead populations are not clear in some channels, be liberal with your gates, 
                  as only events falling in the intersection of all gates will get used for normalization."), 
        fluidPage(
            fluidRow(align="center",
                     # previous sample button
                     div(style="display:inline-block; vertical-align:center",
                         bsButton(inputId="prevSmpl", 
                                  label=NULL,
                                  icon=icon("chevron-left"),
                                  style="default",
                                  size="extra-small")),
                     # sample selection
                     div(style="display:inline-block",  
                         selectInput(
                             inputId="select_sample",
                             label=NULL,
                             choices=samples,
                             selected=samples[[selected]],
                             width="240px")),
                     # next sample button
                     div(style="display:inline-block; vertical-align:center",
                         bsButton(inputId="nextSmpl", 
                                  label=NULL,
                                  icon=icon("chevron-right"),
                                  style="default",
                                  size="extra-small")),
                     div(style="display:inline-block",    
                         h6(strong("Cofactor:"))),
                     # transformation cofactor input
                     div(style="display:inline-block",    
                         numericInput(inputId="input_cfGating", 
                                      label=NULL, 
                                      value=5, 
                                      min=1, 
                                      width="60px")),
                     div(style="display:inline-block",    
                         h6(strong("Number of events:"))),
                     # number of events input
                     div(style="display:inline-block",    
                         numericInput(inputId="input_nGating", 
                                      label=NULL, 
                                      value=20e3, 
                                      min=1e4, 
                                      width="90px")),
                     div(style=inlineTop,
                         bsButton(inputId="button_viewGating", 
                                  label=NULL,
                                  icon=icon("share"),
                                  style="success")),
                     bsTooltip(id="button_viewGating", 
                         title="View",
                         placement="right", 
                         trigger="hover"),
                     div(style=inlineTop,
                         bsButton(inputId="gateBeads", 
                                  label=NULL,
                                  icon=icon("object-ungroup"),
                                  style="primary")),
                     bsTooltip(id="gateBeads", 
                               title="Gate",
                               placement="right", 
                               trigger="hover")
            ),
            # beads vs. dna scatters
            fluidRow(
                align="center",
                plotOutput(outputId="plot_beadScatter1", width="20%",
                           brush=brushOpts(id="beadGate1", resetOnNew=TRUE)),
                plotOutput(outputId="plot_beadScatter2", width="20%",
                           brush=brushOpts(id="beadGate2", resetOnNew=TRUE)),
                plotOutput(outputId="plot_beadScatter3", width="20%",
                           brush=brushOpts(id="beadGate3", resetOnNew=TRUE)),
                plotOutput(outputId="plot_beadScatter4", width="20%",
                           brush=brushOpts(id="beadGate4", resetOnNew=TRUE)),
                plotOutput(outputId="plot_beadScatter5", width="20%",
                           brush=brushOpts(id="beadGate5", resetOnNew=TRUE))
            )
        )
    )
}

# ------------------------------------------------------------------------------
# smoothed beads panel 
# ------------------------------------------------------------------------------
box_smoothedBeads <- 
    shinydashboard::box(
        width=12, 
        collapsible=TRUE, 
        plotOutput("plot_smoothedBeads"))