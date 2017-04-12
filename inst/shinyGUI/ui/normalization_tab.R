inline    <- "display:inline-block;"
inlineTop <- "display:inline-block; vertical-align:top"

# ==============================================================================
# NORMALIZATION 
# ------------------------------------------------------------------------------

collapseBox <- "shinyjs.collapse = function(boxId) {
    $('#'+boxId).closest('.box').find('[data-widget=collapse]').click();}"

normalization_tab <- fluidPage(
    shinyjs::useShinyjs(),
    shinyjs::extendShinyjs(text=collapseBox),
    # ----- styling -----
    tags$head(tags$style(HTML(".small-box{height:96px; 
                                          margin-bottom:0px}"))),
    tags$head(tags$style("#dwnld_normResults{height:34px; 
                                             width:100%; 
                                             margin-top:10px; 
                                             border-width:1.25px; 
                                             border-color:green; 
                                             background-color:honeydew}")),
    tags$head(tags$style("#button_customBeads{height:34px}")),
    tags$head(tags$style("#button_viewGating{height:34px}")),
    tags$head(tags$style("#plot_beadScatter1{float:left}")),
    tags$head(tags$style("#plot_beadScatter2{float:left}")),
    tags$head(tags$style("#plot_beadScatter3{float:left}")),
    tags$head(tags$style("#plot_beadScatter4{float:left}")),
    tags$head(tags$style("#plot_beadScatter5{float:left}")),
    tags$head(tags$style("#plot_smoothedBeads{height:100vh !important;}")),
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

# box 4: gating yield and download results -------------------------------------
box4 <- shinydashboard::box(
    id="box_4",
    width=12, 
    collapsible=TRUE,
    #solidHeader=TRUE,
    valueBox(
        value=textOutput("text_gatingYield"),
        subtitle=textOutput("text_smplNm"),
        icon=icon("percent"),
        color="green", 
        width=NULL),
    downloadButton(
        outputId="dwnld_normResults", 
        label="Download results", 
        style=dwnld_2
    )
)

# ------------------------------------------------------------------------------
# bead gating panel  
# ------------------------------------------------------------------------------
box_beadGating <- function(samples) {
    shinydashboard::box(
        width=12, 
        collapsible=TRUE,
        footer=em(strong("Brush to gate beads."), 
            "Ideally, bead-bead and cell-bead doublets should be excluded.
            However, if bead populations are not clear in some channels, be liberal with your gates, 
            as only events falling in the intersection of all gates will get used for normalization."), 
        fluidPage(
            fluidRow(align="center",
                div(style="display:inline-block; vertical-align:top; height:34px",
                    actionButton(
                        inputId="prevSmpl", 
                        label=NULL,
                        icon=icon("chevron-left"))),
                div(style="display:inline-block; vertical-align:top",  
                    selectInput(
                        inputId="select_sample",
                        label=NULL,
                        choices=samples,
                        selected=samples[[1]],
                        width="240px"
                    )),
                div(style="display:inline-block; vertical-align:top; height:34px",
                    actionButton(
                        inputId="nextSmpl", 
                        label=NULL,
                        icon=icon("chevron-right"))),
                div(style=inline,    h6(strong("Cofactor:"))),
                div(style=inline,    numericInput("input_cfGating", NULL, value=5, min=1, width="60px")),
                div(style=inline,    h6(strong("Number of events:"))),
                div(style=inline,    numericInput("input_nGating", NULL, value=25e3, min=1e4, width="120px")),
                div(style=inlineTop, actionButton("button_viewGating", NULL, icon("eye")))
            ),
            # ----- beads vs. dna scatters -----
            fluidRow(align="center",
                plotOutput("plot_beadScatter1", brush="gate1", width="20%"),
                plotOutput("plot_beadScatter2", brush="gate2", width="20%"),
                plotOutput("plot_beadScatter3", brush="gate3", width="20%"),
                plotOutput("plot_beadScatter4", brush="gate4", width="20%"),
                plotOutput("plot_beadScatter5", brush="gate5", width="20%")
            )
            # -----
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