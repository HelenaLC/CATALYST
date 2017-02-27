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
                               textOutput("text_compCytof_1"),
                               verbatimTextOutput("text_compCytof_2"),
                               uiOutput("compensation_sidebar_1"),
                               uiOutput("compensation_sidebar_2"))))

# ------------------------------------------------------------------------------

compensation_sidebar_1 <- tagList(
    checkboxInput("use_mean", "Use means"),
    checkboxInput("est_trim", "Use trim estimate", value=TRUE),
    checkboxInput("enter_trim", "Enter trim value"),
    uiOutput("enter_trim"),
    actionButton("compensate", "Compensate", style=button_style))

# ------------------------------------------------------------------------------

compensation_sidebar_2 <- tagList(
    hr(style="border-color:black"),
    tags$style(dwnld_1), tags$style(dwnld_2),
    downloadButton("dwnld_1", "Compensated beads", class="dwnld_1"),
    downloadButton("dwnld_2", "Compensated cells", class="dwnld_1"),
    downloadButton("dwnld_3", "Compensation matrix",  class="dwnld_2"))

