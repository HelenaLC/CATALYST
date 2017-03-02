# ------------------------------------------------------------------------------
# compensation tab 
# ------------------------------------------------------------------------------

compensation_tab <- fluidPage(
    tags$head(tags$style("#text_compCytof_1 {color:firebrick}")),
    tags$head(tags$style("#text_compCytof_2 {font-size:x-small}")),
    sidebarLayout(position="left",
                  mainPanel(width=8,
                            tabBox(width=12,
                                   tabPanel(icon("info-circle"),         uiOutput("compensation_guide")),
                                   tabPanel("Trim value estimation",     uiOutput("panel_estTrims")),
                                   tabPanel("Spillover matrix",          uiOutput("panel_plotSpillmat")),
                                   tabPanel("Before vs. after scatters", uiOutput("panel_scatters")),
                                   tabPanel("Summary plot",              uiOutput("panel_plotScatter")))),
                  sidebarPanel(width=4,
                               fileInput("fcs2", "Upload FCS", accept=".fcs"),
                               uiOutput("compensation_sidebar_1"),
                               uiOutput("compensation_sidebar_2"))))

# ------------------------------------------------------------------------------

compensation_sidebar_1 <- tagList(
    checkboxInput("box_estSM", "Estimate spill from single-stained controls", value=TRUE),
    checkboxInput("box_upldSM", "Upload spillover matrix (CSV)"),
    uiOutput("input_upldSM"),
    uiOutput("text_compCytof"), 
    checkboxInput("box_estTrim", "Estimate trim value", value=TRUE),
    checkboxInput("box_enterTrim", "Enter trim value"),
    checkboxInput("box_useMedians", "Use medians"),
    uiOutput("text_enterTrim"),
    div(uiOutput("input_enterTrim"),  style="display:inline-block; width:25%"),
    div(uiOutput("button_enterTrim"), style="display:inline-block; width:25%; vertical-align:top"))

# ------------------------------------------------------------------------------

compensation_sidebar_2 <- tagList(
    hr(style="border-color:black"),
    tags$style(dwnld_1), tags$style(dwnld_2),
    downloadButton("dwnld_comped_1", "Compensated beads", class="dwnld_1"),
    downloadButton("dwnld_comped_2", "Compensated cells", class="dwnld_1"),
    downloadButton("dwnld_spillMat", "Spillover matrix",  class="dwnld_2"))

