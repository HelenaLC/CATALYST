# ------------------------------------------------------------------------------
# compensation tab 
# ------------------------------------------------------------------------------

compensation_tab <- fluidPage(
    tags$head(tags$style("#text_compCytof_1 {color:firebrick}")),
    tags$head(tags$style("#text_compCytof_2 {font-size:x-small}")),
    sidebarLayout(position="left",
        mainPanel(
            width=8,
            tabBox(
                width=12,
                tabPanel(
                    icon("info-circle"), 
                    uiOutput("compensation_guide")),
                tabPanel(
                    "Trim value estimation", 
                    uiOutput("panel_estTrim")),
                tabPanel(
                    "Spillover matrix", 
                    uiOutput("panel_plotSpillmat")),
                tabPanel(
                    "Before vs. after scatters", 
                    uiOutput("panel_scatters")))),
        sidebarPanel(
            width=4,
            fileInput("fcs2", "Upload FCS", accept=".fcs"),
            uiOutput("compensation_sidebar_1"),
            uiOutput("compensation_sidebar_2"))))

# ------------------------------------------------------------------------------

compensation_sidebar_1 <- tagList(
    checkboxInput("box_estSM", 
        "Estimate spill from single-stained controls", value=TRUE),
    checkboxInput("box_upldSM", 
        "Upload spillover matrix (CSV)"),
    uiOutput("input_upldSM"),
    uiOutput("text_compCytof"), 
    h5(strong("Enter trim value for estimation of spill values")),
    div(numericInput("input_enterTrim", NULL,
        value=NULL, min=.01, max=.5, step=.01),  
        style="display:inline-block; width:25%"),
    div(actionButton("button_enterTrim", "Confirm", style=ylw_button), 
        style="display:inline-block; width:25%; vertical-align:top"),
    helpText("Note that a trim value of 0.5 is equal to using medians."))

# ------------------------------------------------------------------------------

compensation_sidebar_2 <- tagList(
    hr(style="border-color:black"),
    tags$style(dwnld_1), tags$style(dwnld_2),
    downloadButton("dwnld_comped_1", "Compensated beads", class="dwnld_1"),
    downloadButton("dwnld_comped_2", "Compensated cells", class="dwnld_1"),
    downloadButton("dwnld_spillMat", "Spillover matrix",  class="dwnld_2"))

