# ==================================================================================================
# DEBARCODING
# --------------------------------------------------------------------------------------------------
debarcoding_guide <- fluidPage(
    tags$style(HTML("li{list-style-type:decimal; padding-bottom:5px; padding-top:5px}")),
    h1("Single-cell deconvolution", style="color:firebrick"),
    
    p("CATALYST provides an implementation of the single-cell deconvolution algorithm described in", 
      a("Zunder et al., 2015:", em("Palladium-based mass tag cell barcoding with a doublet-filtering 
                                   scheme and single-cell deconvolution algorithm."), 
        href="http://www.nature.com/nprot/journal/v10/n2/full/nprot.2015.020.html")),
    
    p("In brief, events are assigned to a sample when 
      i) their positive and negative barcode populations are separated by a distance 
      larger than a threshold value and ii) the combination of their positive barcode channels 
      appears in the provided barcoding scheme. CATALYST optionally allows for", 
      strong("population-specific"), "separation thresholds, as well as", 
      strong("automated estimation"), "and", strong("individual adjustment"), "of these thresholds"),
    
    box(width=12, status="primary", style="background-color:aliceblue", tags$ul(
        tags$li(strong("Upload FCS"), "- the measurement data you wish to debarcode"),
        tags$li(strong("Upload barcoding scheme"), "- a binary matrix 
                with barcode masses as column and sample IDs as row names
                that specifies which channels are to be positive for each sample", 
                HTML("<u>or</u>"), strong("Select barcode channels"), 
                "in the case of single-stained samples",
                div(style="display:inline; color:firebrick", "- these may be used for compensation.")),
        tags$li(strong("Assign preliminary IDs", style="color:steelblue")),
        tags$li(strong("Estimate separation cutoffs"), 
                "(this will be done automatically) or", 
                strong("Enter global separation cutoffs"), 
                "and select a", strong("Mahalanobis distance threshold")), 
        tags$li("Review", strong("Yield, Event and Mahal plot", style="color:dimgrey"), 
                "and", strong("Adjust population-specific cutoffs"), 
                "to modulate deconvolution stringency and cell yields"),
        tags$li(strong("Apply cutoffs", style="color:steelblue"), 
                "to finalize event assignments"),
        tags$li(strong("Output FCS files", style="color:green"), 
                "to obtain population-wise FCS files and", 
                strong("Get plots", style="color:green"), 
                "to download diagnostic plots"))),
              
    h4("Assignment of preliminary barcode IDs", style="color:darkorange"),
    p("The debarcoding step commences by preliminarily assigning each event to a sample.
      In the simplest case of single-stained data, the most intense channel is considered positive 
      and its respective mass assigned as ID. Elsewise, provided a barcoding scheme with a", 
      em("coherent number k of positive channels", style="font-weight:bold"), 
      "for all samples, the k highest channels are considered positive and n−k channels negative.
      Given a", em("non-constant number of 1’s", style="font-weight:bold"), "between samples, 
      the highest separation between consecutive barcodes is looked at. 
      Based on channels considered positive and negative (1 or 0, respectively), 
      each event is assigned a binary code that,", HTML("<u>if</u>"), 
      "matched with a code in the barcoding scheme, dictates which ID it is given.
      Events whose positive barcodes are still very low or whose binary pattern of 
      positive and negative barcodes doesn’t occur in the barcoding scheme 
      will be given ID 0 for “unassigned”."), 

    h4("Application of deconvolution parameters", style="color:darkorange"),
    
    p("The", strong("Yield plot"), "tab displays a histogram of events separated by a given distance, 
      as well as yields upon debarcoding as a function of separation cutoffs. 
      Selecting 0 will render a summary plot of all barcodes. Complementary, the right-hand table 
      summarizes current cutoff choices, event assignments and yields.")
)

# ==================================================================================================
# COMPENSATION
# --------------------------------------------------------------------------------------------------

compensation_guide <- fluidPage(
   
    h1("Bead-based compensation", style="color:firebrick"),
    
    p("CATALYST performs compensation via a two-step approach comprising: 
      i) identification of single positive populations; and, ii) estimation of a spillover matrix
      from the populations identified, followed by compensation via multiplication of 
      measurement intensities by its inverse, the compensation matrix."),

    
    box(width=12, status="primary", style="background-color:aliceblue", tags$ul(
        tags$li(strong("Upload FCS"), "- the measurement data you wish to compensate"),
        tags$li(strong("Upload spillover matrix (CSV)"),
                "- a spillover matrix with column and row names", 
                HTML("<u>or</u>"), strong("Estimate spill from single-stained controls"),
                div(style="display:inline; color:firebrick", 
                    "- this requires having selected barcode channels during debarcoding!")),
        tags$li("Review", strong("Before vs. after scatters", style="color:dimgrey"), 
                "to check current compensation and", 
                strong("Adjust", style="color:orange"), 
                "spill values to align positive and negative populations"),
        tags$li("Download", strong("Compensated beads, cells and spillover matrix", style="color:green")))))






