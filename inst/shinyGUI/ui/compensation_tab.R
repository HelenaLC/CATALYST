# ==============================================================================
# compensation tab
# ==============================================================================

compensationTab <- fluidPage(
    tags$style("#plotSpillmat{height:100vh !important; width:100%}"),
    sidebarLayout(
        position="left",
        mainPanel(
            width=9,
            tabBox(
                width=12,
                tabPanel(
                    title=icon("info-circle"), 
                    uiOutput("compensation_guide")),
                tabPanel(
                    title="Spillover matrix", 
                    plotOutput("plotSpillmat")),
                tabPanel(
                    title="Before vs. after scatters", 
                    uiOutput("panel_scatters")))),
        sidebarPanel(
            width=3,
            fileInput(
                inputId="fcsComp", 
                label="Upload FCS", 
                multiple=TRUE,
                accept=".fcs"),
            uiOutput("compensationSidebar1"),
            uiOutput("compensationSidebar2"))))

# ------------------------------------------------------------------------------
# 1st sidebar
# ------------------------------------------------------------------------------

compensationSidebar1 <- tagList(
    checkboxInput(
        inputId="box_upldSm", 
        label="Upload spillover matrix (CSV)"),
    uiOutput("inputSm"),
    checkboxInput(
        inputId="box_estSm", 
        label="Estimate spill from single-stained controls"))

# ------------------------------------------------------------------------------
# 2nd sidebar
# ------------------------------------------------------------------------------

compensationSidebar2 <- tagList(
    hr(style="border-color:black"),
    checkboxInput(
        inputId="box_setToZero", 
        label="Should negative values be set to zero?"),
    tags$style(type="text/css", "#dwnld_comped {
        display:inline-block; color:white; width:49%; float:left}"),
    tags$style(type="text/css", "#dwnld_spillMat {
        display:inline-block; color:white; width:49%; float:right}"),
    downloadButton(
        outputId="dwnld_comped", 
        label="Compensated data",
        class="btn-success"), 
    downloadButton(
        outputId="dwnld_spillMat", 
        label="Spillover matrix", 
        class="btn-success"))
