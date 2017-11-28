# ==============================================================================
# user guide tab
# ==============================================================================

guidesTab <- fluidPage( 
    tags$style(type="text/css", "a{color:steelblue; font-weight:bold}"),
    sidebarLayout(
        mainPanel(
            width=9, 
            tabBox(
                width=12, 
                tabPanel(
                    title="Concatenation", 
                    uiOutput("concatenation_guide")),
                tabPanel(
                    title="Normalization", 
                    uiOutput("normalization_guide")),

                tabPanel(
                    title="Compensation", 
                    uiOutput("compensation_guide")),
                tabPanel(
                    title="Debarcoding", 
                    uiOutput("debarcoding_guide")))),
        sidebarPanel(
            width=3,
            "We designed CATALYST (Cytometry dATa anALYSis Tools) 
            to provide tools for preprocessing of cytometry data, including 
            i) concatenation, ii) normalization using bead standards, 
            iii) single-cell deconvolution, and iv) bead-based compensation.",
            strong(style="color:firebrick", "The preprint to this project 
                is accesible on", strong("bioRxiv"), "at",
                a("Chevrier et al., 2017:", em("Channel crosstalk correction 
                in suspension and imaging mass cytometry."), 
                href="https://www.biorxiv.org/content/early/2017/09/07/185744")))))

# ------------------------------------------------------------------------------
# concatenation
# ------------------------------------------------------------------------------
concatenation_guide <- fluidPage(
    h1("FCS file concatenation", 
        style="color:firebrick"),
    shinydashboard::box(
        width=12, 
        status="primary", 
        style="background-color:aliceblue", 
        tags$ul(
            tags$li(strong("Upload FCS", style="color:orange"), 
                "- the measurement data you wish to concatenate.",
                div(style="color:firebrick",
                    "(note that all acquisitions must contain 
                    identical measurement parameters!)")),
            tags$li("Specify the", strong("Output file name")),
            tags$li("Specify whether to", strong("Add file number channel"), 
                "and whether to", strong("Order by time"), "of acquisition"), 
            tags$li("Optional", strong("File editing", style="color:orange"), 
                "allows modification of the common 
                parameter descriptions of the output FCS file"),
            tags$li(strong("Go to normalization", style="color:dimgrey"), 
               "if you want to directly propagate the data to the next
                preprocessing step"),
            tags$li(strong("Merges files", style="color:green"),
                "to download a single FCS file of the concatenated data"))))

# ------------------------------------------------------------------------------
# normalization
# ------------------------------------------------------------------------------
normalization_guide <- fluidPage(
    h1("Normalization using bead standards", 
        style="color:firebrick"),
    p(a("Finck et al., 2013:", 
        em("Normalization of mass cytometry data with bead standards."), 
        href="http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.22271/abstract")),
    shinydashboard::box(
        width=12, 
        status="primary", 
        style="background-color:aliceblue", 
        tags$ul(
            tags$li(strong("Upload FCS", style="color:orange"), 
                "- this may be a single or multiple FCS file(s) 
                of the measurement data you wish to normalize", 
                div(style="color:dimgrey", 
                    "(Note that every sample will require gating! 
                    If you have many samples, it is recommended to 
                    concatenate them prior to normalization.)")),
            tags$li(strong("Normalization baseline", style="color:orange"), 
                "- samples may be normalized to either the median level 
                of the input data, or to that of previously exported beads"), 
            tags$li(strong("Bead selection", style="color:orange"), 
                "- select which mass channels to use as standards."),
            tags$li(strong("Bead gating", style="color:dimgrey"), 
                "- once beads have been selected, bead vs. DNA scatters will 
                be rendered, and you may brush points to gate bead populations. 
                When all gates have been drawn, finalize them by pressing",
                strong("Gate", style="color:steelblue")),
            tags$li(strong("Bead removal", style="color:dimgrey"), "- when", 
                strong("Should beads be removed?"), "is checked, bead vs. bead 
                scatters color coded by mahalanobis distance from identified 
                beads will be displayed. Select a cutoff below which events 
                should be remove by dragging the slider and",
                strong("Apply", style="color:steelblue"), "it"),
            tags$li(strong("Propagating / downloading data"),
                "- the right most box will display the percentage of events 
                that were identified as beads and the fraction of samples 
                that have been gated. You may directly", 
                strong("Go to compensation", style="color:dimgrey"), 
                "or download the", 
                strong("Normalized data", style="color:green")))))

# ------------------------------------------------------------------------------
# debarcoding
# ------------------------------------------------------------------------------
debarcoding_guide <- fluidPage(
    tags$style(HTML("li{
        list-style-type:square; 
        padding-bottom:5px; 
        padding-top:5px}")),
    h1("Single-cell deconvolution", 
        style="color:firebrick"),
    p(a("Zunder et al., 2015:", em("Palladium-based mass tag cell barcoding 
        with a doublet-filtering scheme and single-cell deconvolution algorithm."), 
        href="http://www.nature.com/nprot/journal/v10/n2/full/nprot.2015.020.html")),
    shinydashboard::box(
        width=12, 
        status="primary", 
        style="background-color:aliceblue", 
        tags$ul(
            tags$li(strong("Upload FCS"), "- the measurement data you wish to debarcode"),
            tags$li(strong("Upload barcoding scheme (CSV)"), "- a binary matrix 
                with barcode masses as column and sample IDs as row names
                that specifies which channels are to be positive for each sample"),
            tags$li(strong("Debarcode", style="color:orange"), 
                "to assign a preliminary ID to each event"),
            tags$li(strong("Estimate separation cutoffs"), 
                "(this will be done automatically) or", 
                strong("Enter global separation cutoff")),
            tags$li("Review", strong("Yield plot", style="color:dimgrey"),
                "and", strong("Adjust population-specific cutoffs"), 
                "to modulate deconvolution stringency and cell yields"),
            tags$li("Review", strong("Mahal plot", style="color:dimgrey"), 
                "and apply a", strong("Mahalanobis distance threshold")),
            tags$li("Download population-wise", 
                strong("FCS files", style="color:green"), "and diagnostic",
                strong("Plots", style="color:green")))))

# ------------------------------------------------------------------------------
# compensation
# ------------------------------------------------------------------------------
compensation_guide <- fluidPage(
    h1("Bead-based compensation", 
        style="color:firebrick"),
    p(a("Chevrier et al., 2017:", em("Channel crosstalk correction 
        in suspension and imaging mass cytometry."), 
        href="https://www.biorxiv.org/content/early/2017/09/07/185744")),
    shinydashboard::box(
        width=12, 
        status="primary", 
        style="background-color:aliceblue", 
        tags$ul(
            tags$li(strong("o Use pre-acquired spillover matrix"),
                "- a spillover matrix with column and row names", 
                HTML("<u>or</u>"), strong("o Estimate spill from controls")),
            tags$li("To estimate spill from controls,",
                strong("Upload single-stains (FCS)"), ",",
                strong("Select single-positive channels"), "and",
                strong("Deconvolute single-stains", style="color:orange")),
            tags$li(strong("Select method"), 
                "to use to compensate the data:",
                strong("o NNLS compensation"), HTML("<u>or</u>"),
                "classical", strong("o Flow compensation")),
            tags$li(strong("Upload multiplexed data (FCS)"), 
                "- the measurement data you wish to compensate"),
            tags$li("Review", 
                strong("Before vs. after scatters", style="color:dimgrey"), 
                "to check the current compensation and", 
                strong("Adjust", style="color:orange"), 
                "spill values to align positive and negative populations"),
            tags$li("Download the", 
                strong("Compensated data", style="color:green"), "and/or",
                strong("Spillover matrix", style="color:green")),
            tags$li(strong("Go to debarcoding", style="color:dimgrey"), 
                "to directly propagate the compensated data"))))