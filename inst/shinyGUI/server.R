# set maximum web request size to 5 GB
options(shiny.maxRequestSize=5000*1024^2)

shinyServer(function(input, output, session) {
    
    output$concatenation_guide <- renderUI(concatenation_guide)
    output$normalization_guide <- renderUI(normalization_guide)
    output$compensation_guide  <- renderUI(compensation_guide)
    output$debarcoding_guide   <- renderUI(debarcoding_guide)
    
    source("server-concatenation.R", local=TRUE)
    source("server-normalization.R", local=TRUE)
    source("server-compensation.R",  local=TRUE)
    source("server-debarcoding.R",   local=TRUE)
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
        mhlCutoffComp = 30   # default Mahalanobis distance cutoff
    )
    if (!exists("custom_isotope_list")) {
        showNotification(closeButton=FALSE,
            h4(strong("Using default isotope list.")))
        vals$isotope_list <- CATALYST::isotope_list
    } else {
        # validate custom isotope list:
        # 1. is list
        # 2. all elements are named
        # 3. no elements are empty
        # 4. all elements are numeric
        # 5. no elements contain duplicates
        l <- custom_isotope_list
        check1 <- !is.list(l)
        check2 <- length(names(l)) != length(l)
        check3 <- any(sapply(l, length) == 0)
        check4 <- !all(sapply(l, is.numeric))
        check5 <- any(!equals(sapply(l, length), 
            sapply(l, function(x) length(unique(x)))))
        if (any(c(check1, check2, check3, check4, check5))) {
            showNotification(duration=NULL,
                "Specified custom isotope list is invalid.
                Please make sure `custom_isotope_list` is a fully
                named list with unique, positive numeric elements.")
            showNotification(closeButton=FALSE,
                h4(strong("Using default isotope list.")))
            vals$isotope_list <- CATALYST::isotope_list
            return()
        }
        showNotification(closeButton=FALSE,
            h4(strong("Using custom isotope list.")))
        vals$isotope_list <- custom_isotope_list
    }
})