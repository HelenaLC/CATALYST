# ==================================================================================================
# DEBARCODING
# --------------------------------------------------------------------------------------------------
debarcoding_guide <- fluidPage(
    tags$style(HTML("li{
        list-style-type:decimal; 
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
                that specifies which channels are to be positive for each sample", 
                HTML("<u>or</u>"), strong("Select single-positive channels"),
                "in the case of single antibody stained samples",
                div(style="display:inline; color:firebrick", 
                    "- these may be used to estimate a spillover matrix")),
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

# ==================================================================================================
# COMPENSATION
# --------------------------------------------------------------------------------------------------
compensation_guide <- fluidPage(
    h1("Bead-based compensation", 
        style="color:firebrick"),
    shinydashboard::box(
        width=12, 
        status="primary", 
        style="background-color:aliceblue", 
        tags$ul(
            tags$li(strong("Upload FCS"), 
                "- the measurement data you wish to compensate"),
            tags$li(strong("Upload spillover matrix (CSV)"),
                "- a spillover matrix with column and row names", 
                HTML("<u>or</u>"), 
                strong("Estimate spill from single-stained controls"),
                div(style="display:inline; color:firebrick", "- this requires 
                    having selected single-positive channels during debarcoding!")),
            tags$li("Review", strong("Before vs. after scatters", style="color:dimgrey"), 
                "to check the current compensation and", 
                strong("Adjust", style="color:orange"), 
                "spill values to align positive and negative populations"),
            tags$li("Optionally, set negative values to zero, and download the", 
                strong("Compensated data", style="color:green")))))






