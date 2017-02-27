# ------------------------------------------------------------------------------
# debarcoding tab
# ------------------------------------------------------------------------------

debarcoding_tab <- fluidPage(      
    tags$style(type="text/css", "a{color:steelblue; 
                                   font-weight:bold}"),
    sidebarLayout(
        mainPanel(width=8, 
                  fluidPage(
                      fluidRow(
                          tabBox(width=12,
                              tabPanel(icon("info-circle"), uiOutput("debarcoding_guide")),
                              tabPanel("Yield plot", uiOutput("yp_panel")),
                              tabPanel("Event plot", uiOutput("ep_panel")),
                              tabPanel("Mahal plot", uiOutput("mhl_panel")))))),
        sidebarPanel(width=4,
                     fileInput("fcs", "Upload FCS", accept=".fcs", width="100%"),
                     uiOutput("debarcoding_sidebar_1"),
                     uiOutput("debarcoding_sidebar_2"))))
                                 
                  

# ------------------------------------------------------------------------------

debarcoding_sidebar_1 <- tagList(
    checkboxInput("box_csv", "Upload barcoding scheme (CSV)", value=TRUE),   
    uiOutput("file_csv"),
    checkboxInput("box_bcChs", "Select barcode channels"), 
    uiOutput("select_bcChs"),
    tags$style(type="text/css", "a{border-width:1.25px}"),
    actionButton("button_assignPrelim", "Assign preliminary IDs", style=button_style))
                                         
# ------------------------------------------------------------------------------

debarcoding_sidebar_2 <- tagList(
    hr(style="border-color:black"),
    # cutoff estimates
    checkboxInput("box_estCutoffs", "Estimate separation cutoffs", value=TRUE),
    # population-specific cutoffs
    checkboxInput("box_indivCutoffs", "Adjust population-specific cutoffs"),
    div(uiOutput("select_adjustCutoff"), style="display:inline-block; width:25%"),
    div(uiOutput("input_adjustCutoff"),  style="display:inline-block; width:25%; vertical-align:top"),
    div(uiOutput("button_adjustCutoff"), style="display:inline-block; width:25%; vertical-align:top"),
    # global sparation cutoff
    checkboxInput("box_globalCutoff", "Enter global separation cutoff"),
    div(uiOutput("input_globalCutoff"),  style="display:inline-block; width:25%"),
    div(uiOutput("button_globalCutoff"), style="display:inline-block; width:25%; vertical-align:top"),
    # mhl cutoff slider
    sliderInput("input_mhlCutoff", "Mahalanobis distance threshold",
                width="100%", min=0, max=100, value=30),
    # "Apply cutoffs" button
    actionButton("button_applyCutoffs", "Apply cutoffs", style=button_style),
    hr(style="border-color:black"),
    # download buttons
    tags$style(dwnld_1), tags$style(dwnld_2),
    downloadButton("dwnld_fcs", "Output FCS files", class="dwnld_1"), 
    downloadButton("dwnld_yep", "Get plots", class="dwnld_2"))




            