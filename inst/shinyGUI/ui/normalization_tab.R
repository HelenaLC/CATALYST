# ==================================================================================================
# NORMALIZATION 
# --------------------------------------------------------------------------------------------------

normalization_tab <- fluidPage(
    tags$head(tags$style("#text_compCytof_1 {color:firebrick}")),
    tags$head(tags$style("#text_compCytof_2 {font-size:x-small}")),
    sidebarLayout(mainPanel(width=8,
                            tabBox(width=12,
                                   tabPanel(icon("info-circle"),   uiOutput("normalization_guide")),
                                   tabPanel("Bead identification", uiOutput("panel_getBeads")))),
                  sidebarPanel(width=4,
                               fileInput("normalization_fcs", "Upload FCS", accept=".fcs"),
                               uiOutput("normalization_sidebar_1"),
                               uiOutput("normalization_sidebar_2"))))

# --------------------------------------------------------------------------------------------------

normalization_sidebar_1 <- tagList(
    selectInput("select_beadMasses", NULL, choices=c("DVS Beads (140, 151, 153, 165, 175)",
                                                     "Beta Beads (139, 141, 159, 169, 175",
                                                     "Custom")),
    uiOutput("input_beadMasses"),
    uiOutput("text_compCytof"), 
    checkboxInput("box_autoBeads", "Identify beads automatically", value=TRUE),
    checkboxInput("box_gateBeads", "Gate beads manually"))

# --------------------------------------------------------------------------------------------------

normalization_sidebar_2 <- tagList(
    hr(style="border-color:black"),
    tags$style(dwnld_1), 
    downloadButton("dwnld_normedBeads", "Normalized data", class="dwnld_1"))

