# set maximum web request size to 5 GB
options(shiny.maxRequestSize=5000*1024^2)

shinyServer(function(input, output, session) {
    
    #source("helpers.R")
    source("guides.R")
    output$debarcoding_guide  <- renderUI(debarcoding_guide)
    output$compensation_guide <- renderUI(compensation_guide)

    source("server-concatenation.R", local=TRUE)
    source("server-normalization.R", local=TRUE)
    source("server-debarcoding.R",   local=TRUE)
    source("server-compensation.R",  local=TRUE)
    source("server-FCSchecking.R",   local=TRUE)

# ------------------------------------------------------------------------------
    
    vals <- reactiveValues(
        # logicals indicating whether data should be propagated
        keepDataConcat = FALSE,
        keepDataNorm = FALSE,
        keepDataComp = FALSE,
        # debarcoding
        debaKeyIsValid = FALSE, # set to TRUE if input CSV passes validity check
        dbFrame1Deba = NULL,    # preliminary dbFrame
        dbFrame2Deba = NULL,    # dbFrame with deconvolution parameters applied
        mhlCutoffDeba = 30,     # default Mahalanobis distance cutoff
        # compensation
        dbFrame1Comp = NULL, # preliminary dbFrame
        dbFrame2Comp = NULL, # dbFrame with deconvolution parameters applied
        mhlCutoffComp = 30   # default Mahalanobis distance cutoff
    )
})