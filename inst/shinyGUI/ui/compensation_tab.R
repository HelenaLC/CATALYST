# ------------------------------------------------------------------------------
# compensation tab 
# ------------------------------------------------------------------------------

compensation_tab <- fluidPage(
    tags$head(tags$style("#text_compCytof_1 {color:firebrick}")),
    tags$head(tags$style("#text_compCytof_2 {font-size:x-small}")),
    sidebarLayout(position="left",
                  mainPanel(width=8,
                            tabBox(width=12,
                                   tabPanel(icon("info-circle"), uiOutput("compensation_guide")),
                                   tabPanel("Estimate trim",    uiOutput("panel_estTrims")),
                                   tabPanel("Spillover matrix", uiOutput("panel_plotSpillmat")),
                                   tabPanel("Summary plot",     uiOutput("panel_plotScatter")),
                                   tabPanel("Before vs. after", uiOutput("panel_scatters")))
                            ),
                  sidebarPanel(width=4,
                               fileInput("fcs2", "Upload FCS", accept=".fcs"),
                               uiOutput("compensation_sidebar_1"),
                               uiOutput("compensation_sidebar_2"))))

# ------------------------------------------------------------------------------

compensation_sidebar_1 <- tagList(
    checkboxInput("box_upldSM", "Upload spillover matrix (CSV)", value=TRUE),
    uiOutput("input_upldSM"),
    checkboxInput("box_estSM",    "Estimate spill from single-stained controls"),
    uiOutput("text_compCytof"), 
    checkboxInput("box_useMedians", "Use medians"),
    checkboxInput("box_estTrim", "Use trim estimate", value=TRUE),
    checkboxInput("box_enterTrim", "Enter trim value"),
    uiOutput("input_enterTrim"))
    #actionButton("button_compensate", "Compensate", style=button_style))

# ------------------------------------------------------------------------------

compensation_sidebar_2 <- tagList(
    hr(style="border-color:black"),
    tags$style(dwnld_1), tags$style(dwnld_2),
    downloadButton("dwnld_comped_1", "Compensated beads", class="dwnld_1"),
    downloadButton("dwnld_comped_2", "Compensated cells", class="dwnld_1"),
    downloadButton("dwnld_spillMat", "Spillover matrix",  class="dwnld_2"))

