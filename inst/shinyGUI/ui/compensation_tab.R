# ==============================================================================
# compensation tab
# ==============================================================================

compensationTab <- fluidPage(
    tags$head(tags$style("#text_compCytof_1 {color:firebrick}")),
    tags$head(tags$style("#text_compCytof_2 {font-size:x-small}")),
    sidebarLayout(position="left",
        mainPanel(
            width=9,
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
            width=3,
            fileInput(inputId="fcsComp", 
                      label="Upload FCS", 
                      multiple=TRUE,
                      accept=".fcs"),
            uiOutput("compensationSidebar1"),
            uiOutput("compensationSidebar2"))))

# ------------------------------------------------------------------------------
# 1st sidebar
# ------------------------------------------------------------------------------

compensationSidebar1 <- tagList(
    checkboxInput(inputId="box_estSm", 
                  label="Estimate spill from single-stained controls"),
    uiOutput("enterTrim"),
    checkboxInput(inputId="box_upldSm", 
                  label="Upload spillover matrix (CSV)"),
    uiOutput("inputSm"),
    uiOutput("text_compCytof"))

enterTrim <- tagList(
    h5(strong("Enter trim value for estimation of spill values")),
    div(style="display:inline-block; width:25%",
        numericInput(inputId="input_enterTrim", 
                     label=NULL,
                     value=NULL, 
                     min=.01, 
                     max=.5, 
                     step=.01)),
    div(style="display:inline-block; width:25%; vertical-align:top",
        shinyBS::bsButton(inputId="button_enterTrim",
                          label=NULL,
                          icon=icon("share"),
                          style="success")),
    shinyBS::bsTooltip(id="button_enterTrim",
                       title="Estimate spill",
                       placement="right",
                       trigger="hover"),
    helpText("Note that a trim value of 0.5 is equal to using medians.")
)

# ------------------------------------------------------------------------------
# 2nd sidebar
# ------------------------------------------------------------------------------

compensationSidebar2 <- tagList(
    hr(style="border-color:black"),
    tags$style(dwnld_1), 
    tags$style(dwnld_2),
    downloadButton(outputId="dwnld_comped_1", 
                   label="Compensated beads", 
                   class="dwnld_1"),
    downloadButton(outputId="dwnld_comped_2", 
                   label="Compensated cells", 
                   class="dwnld_1"),
    downloadButton(outputId="dwnld_spillMat", 
                   label="Spillover matrix",  
                   class="dwnld_2"))
