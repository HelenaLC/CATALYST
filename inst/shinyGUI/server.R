# set maximum web request size to 500 MB
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
    source("helpers.R")
    
# guides 
    source("ui/guides.R")
    output$debarcoding_guide  <- renderUI(debarcoding_guide)
    output$compensation_guide <- renderUI(compensation_guide)

    source("server-norm.R", TRUE)
    source("server-deba.R", TRUE)
    source("server-comp.R", TRUE)
# plot files -----    
    source("ui/debarcoding_plots.R")
    source("ui/compensation_plots.R")

# ------------------------------------------------------------------------------
    
    vals <- reactiveValues(
# debarcoding
        
        dbFrame1 = NULL, # preliminary dbFrame
        dbFrame2 = NULL, # dbFrame with deconvolution parameters applied
        mhlCutoffDeba = 30, # default Mahalanobis distance cutoff
# compensation
        cmp1 = NULL, # compensated flowFrame 1
        cmp2 = NULL) # compensated flowFrame 2
})
