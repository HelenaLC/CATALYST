# ==================================================================================================
# NORMALIZATION 
# --------------------------------------------------------------------------------------------------

normalization_tab <- fluidPage(
    tags$head(tags$style("#text_compCytof_1 {color:firebrick}")),
    tags$head(tags$style("#text_compCytof_2 {font-size:x-small}")),
    sidebarLayout(mainPanel(width=8,
                            tabBox(width=12,
                                   tabPanel(icon("info-circle"),   uiOutput("normalization_guide")),
                                   tabPanel("Bead identification", uiOutput("gating_panel")))),
                  sidebarPanel(width=4,
                               fileInput("fcs_norm", "Upload FCS", accept=".fcs"),
                               uiOutput("normalization_sidebar_1"),
                               uiOutput("normalization_sidebar_2"))))

# --------------------------------------------------------------------------------------------------

normalization_sidebar_1 <- tagList(
    selectInput("select_beads", NULL, 
                choices=c("DVS Beads (140, 151, 153, 165, 175)" = "dvs",
                          "Beta Beads (139, 141, 159, 169, 175)" = "beta",
                          "Custom" = "custom")),
    uiOutput("input_beadMs"),
    uiOutput("text_customBeads"), 
    uiOutput("input_customBeads"), 
    checkboxInput("box_autoBeads", "Identify beads automatically"),
    checkboxInput("box_gateBeads", "Gate beads manually"),
    actionButton("button_gateBeads", "Go", style=button_style))

# --------------------------------------------------------------------------------------------------

normalization_sidebar_2 <- tagList(
    hr(style="border-color:black"),
    tags$style(dwnld_1), 
    downloadButton("dwnld_normedBeads", "Normalized data", class="dwnld_1"))

panel_getBeads <- fluidPage(
    tags$head(tags$style("#plot_bead_identification{height:50vh !important;}")),
    plotOutput("plot_getBeads", width="100%"))
