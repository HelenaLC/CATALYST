# ==============================================================================
# debarcoding tab
# ==============================================================================

debarcodingTab <- fluidPage( 
    tags$style(type="text/css", "a{color:steelblue; font-weight:bold}"),
    shinyjs::useShinyjs(),
    sidebarLayout(
        mainPanel(
            width=9, 
            tabBox(
                width=12, 
                tabPanel(title=icon("info-circle"), 
                         uiOutput("debarcoding_guide")),
                tabPanel(title="Yield plot", 
                         uiOutput("yp_panel")),
                tabPanel(title="Event plot", 
                         uiOutput("ep_panel")),
                tabPanel(title="Mahal plot", 
                         uiOutput("mhl_panel")))),
        sidebarPanel(
            width=3,
            fileInput(inputId="fcsDeba", 
                      label="Upload FCS", 
                      accept=".fcs", 
                      width="100%"),
            uiOutput("debarcodingSidebar1"),
            uiOutput("debarcodingSidebar2"),
            uiOutput("debarcodingSidebar3"))
    )
)

# ------------------------------------------------------------------------------
# 1st sidebar
# ------------------------------------------------------------------------------

debarcodingSidebar1 <- tagList(
    checkboxInput(inputId="box_csv", 
                  label="Upload barcoding scheme (CSV)", 
                  value=TRUE),   
    uiOutput("file_csv"),
    checkboxInput(inputId="box_bcChs", 
                  label="Select single-positive channels"), 
    uiOutput("select_bcChs"),
    shinyBS::bsButton(inputId="button_assignPrelim", 
                      label="Debarcode", 
                      style="primary",
                      block=TRUE))

# ------------------------------------------------------------------------------
# 2nd sidebar
# ------------------------------------------------------------------------------

debarcodingSidebar2 <- tagList(
    hr(style="border-color:black"),
    # automated cutoff estimation
    checkboxInput(inputId="box_estCutoffs", 
                  label="Estimate separation cutoffs", 
                  value=TRUE),
    # population-specific cutoffs
    checkboxInput(inputId="box_adjustCutoff", 
                  label="Adjust population-specific cutoffs"),
    div(style="display:inline-block; width:25%",
        uiOutput("select_adjustCutoff")),
    div(style="display:inline-block; width:25%",
        uiOutput("input_adjustCutoff")),
    div(style="display:inline-block; vertical-align:top",
        uiOutput("button_adjustCutoff")),
    # global sparation cutoff
    checkboxInput(inputId="box_globalCutoff", 
                  label="Enter global separation cutoff"),
    div(style="display:inline-block; width:25%",
        uiOutput("input_globalCutoff")),
    div(style="display:inline-block; vertical-align:top",
        uiOutput("button_globalCutoff")),
    # mahalanobis distance cutoff
    div(style="display:inline-block; vertical-align:middle; width:75%",
        sliderInput(inputId="slider_mhlCutoff", 
                    label="Mahalanobis distance threshold",
                    width="100%", min=10, max=100, value=30)),
    div(style="display:inline-block; vertical-align:middle", 
        shinyBS::bsButton(inputId="button_mhlCutoff",
                          label=NULL,
                          icon=icon("share"),
                          style="warning")),
    shinyBS::bsTooltip(id="button_mhlCutoff",
                       title="Apply",
                       placement="right"),
    hr(style="border-color:black"),
    # checkboxes for output FCS file naming
    checkboxInput(inputId="box_IDsAsNms", 
                  label="Use sample IDs as files names", 
                  value=TRUE),
    checkboxInput(inputId="box_upldNms",  
                  label="Upload naming sheet (CSV)"),
    uiOutput("upldNms"),
    # download buttons
    tags$style(type="text/css", 
               "#dwnld_fcs {display:inline-block; color:white; width:25%}"),
    tags$style(type="text/css", 
               "#dwnld_yep {display:inline-block; color:white; width:25%}"),
    downloadButton(outputId="dwnld_fcs", 
                   label="FCS files",
                   class="btn-success"), 
    downloadButton(outputId="dwnld_yep", 
                   label="Plots", 
                   class="btn-success"))

