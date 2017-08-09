# ------------------------------------------------------------------------------
# check input FCS files
# ------------------------------------------------------------------------------

# get named list of indices of duplicate masses
# in single- and multiplexed data
duplicateMsSS <- reactive(get_duplicate_ms(ssChs()$Mass))
duplicateMsMP <- reactive(get_duplicate_ms(mpChs()$Mass))

# ------------------------------------------------------------------------------
# 1st modal: "Duplicate masses detected" in single-stains
# ------------------------------------------------------------------------------
observe({
    # make sure this pops up only once
    req(is.null(vals$ffControls_msChecked))
    # no selection is required if all masses are unique
    inds <- duplicateMsSS()
    n <- length(inds)
    if (n == 0) {
        vals$ffControls_msChecked <- ffControls()
        return()
    }
    # render modal
    showModal(editMs(id="SS"))
    # render list of masses and selectInputs for metal selection
    output$duplicateMassesSS <- renderUI(msList(n=n, id="SS"))
    output$duplicateMetalsSS <- renderUI(selectMs(
        mets=ssChs()$Metal, inds=inds, id="SS"))
    # render textOutputs of masses
    sapply(seq_len(n), function(i) {
        output[[paste0("massSS", i)]] <- renderText(names(inds)[i])
    })
})

# ------------------------------------------------------------------------------
# 2nd modal: "Metal mismatches detected" b/w Fluidigm & single-stains
# ------------------------------------------------------------------------------
observe({
    # make sure this pops up only once
    req(is.null(vals$ffControls_metsChecked))
    # no checking is required if metals are identical 
    # b/w the Fluidigm default and input single-stains
    if (sum(mismatchingMetsSS()) == 0) {
        vals$ffControls_metsChecked <- ffControls()
        return()  
    }
    inds <- matchingMsSS()
    n <- sum(inds[[1]])
    # render modal
    showModal(editMets(
        id="FLvsSS", 
        header1="Fluidigm", 
        header2="Single-stains"))
    # render list of masses and metals for matching channels
    output$massesFLvsSS  <- renderUI(msList(n=n, id="FLvsSS"))
    output$metals1FLvsSS <- renderUI(metsList(n=n, id1="FLvsSS", id2="ref"))
    output$metals2FLvsSS <- renderUI(metsList(n=n, id1="FLvsSS", id2="dat"))
    # render checkboxInputs for metal selection
    output$boxes1FLvsSS <- renderUI(boxList(n=n, id1="FLvsSS", id2="ref"))
    output$boxes2FLvsSS <- renderUI(boxList(n=n, id1="FLvsSS", id2="dat"))
    # render textOutputs of masses and metals
    sapply(seq_len(n), function(i) {
        id1 <- paste0("massFLvsSS", i)
        id2 <- paste0("metalFLvsSSref", i)
        id3 <- paste0("metalFLvsSSdat", i)
        output[[id1]] <- renderText(flChs()$Mass[inds$fl][i])
        output[[id2]] <- renderText(flChs()$Metal[inds$fl][i])
        output[[id3]] <- renderText(ssChs()$Metal[inds$ss][i])
    })
})

# ------------------------------------------------------------------------------
# 3rd modal: "Duplicate masses detected" in multiplexed data
# ------------------------------------------------------------------------------
observe({
    # make sure this pops up only once
    req(is.null(vals$fsComp_msChecked))
    # no selection is required if all masses are unique
    inds <- duplicateMsMP()
    n <- length(inds)
    if (n == 0) {
        vals$fsComp_msChecked <- fsComp()
        return()
    }
    # render modal
    showModal(editMs(id="MP"))
    # render list of masses and selectInputs for metal selection
    output$duplicateMassesMP <- renderUI(msList(n=n, id="MP"))
    output$duplicateMetalsMP <- renderUI(selectMs(
        mets=mpChs()$Metal, inds=inds, id="MP"))
    # render textOutputs of masses
    sapply(seq_len(n), function(i) {
        output[[paste0("massMP", i)]] <- renderText(names(inds)[i])
    })
})

# ------------------------------------------------------------------------------
# 4th modal: "Metal mismatches detected" b/w single- & multiplexed data
# ------------------------------------------------------------------------------
observe({
    # make sure this pops up only once
    req(is.null(vals$fsComp_metsChecked))
    # no checking is required if metals are identical 
    # b/w the Fluidigm default and input single-stains
    if (sum(mismatchingMetsMP()) == 0) {
        vals$fsComp_metsChecked <- fsComp()
        return()  
    }
    inds <- matchingMsMP()
    n <- sum(inds[[1]])
    # render modal
    showModal(editMets(
        id="SSvsMP", 
        header1="Single-stains", 
        header2="Multiplexed"))
    # render list of masses and metals for matching channels
    output$massesSSvsMP  <- renderUI(msList(n=n, id="SSvsMP"))
    output$metals1SSvsMP <- renderUI(metsList(n=n, id1="SSvsMP", id2="ref"))
    output$metals2SSvsMP <- renderUI(metsList(n=n, id1="SSvsMP", id2="dat"))
    # render checkboxInputs for metal selection
    output$boxes1SSvsMP <- renderUI(boxList(n=n, id1="SSvsMP", id2="ref"))
    output$boxes2SSvsMP <- renderUI(boxList(n=n, id1="SSvsMP", id2="dat"))
    # render textOutputs of masses and metals
    sapply(seq_len(n), function(i) {
        id1 <- paste0("massSSvsMP", i)
        id2 <- paste0("metalSSvsMPref", i)
        id3 <- paste0("metalSSvsMPdat", i)
        output[[id1]] <- renderText(ssChs()$Mass[inds$ss][i])
        output[[id2]] <- renderText(ssChs()$Metal[inds$ss][i])
        output[[id3]] <- renderText(mpChs()$Metal[inds$mp][i])
    })
})

# ------------------------------------------------------------------------------
# 1st modal's bsButton:
# remove duplicate masses from single-stains
observeEvent(input$msCheckedSS, {
    inds <- duplicateMsSS()
    remove <- NULL
    for (i in seq_along(inds)) {
        choices <- ssChs()$Metal[inds[[i]]]
        selected <- input[[paste0("metalSelectionSS", i)]]
        remove <- append(remove, which(inds[[i]])[choices != selected])
    }
    ff <- ffControls()
    vals$ffControls_msChecked <- ff[, -remove]
    removeModal()
})

# ------------------------------------------------------------------------------
# 2nd modal's bsButton:
# update metals of single-stains
observeEvent(input$metsCheckedFLvsSS, {
    inds <- matchingMsSS()
    n <- sum(inds[[1]])
    
    ff <- vals$ffControls_msChecked
    # get parameters
    ssPars <- colnames(ff)
    flPars <- readRDS('data/Fluidigm_default_metals.rds')
    for (i in seq_len(n)) {
        # determine if Fluidigm or single-stain metal has been selected
        selected <- which(c(
            input[[paste0("boxFLvsSSref", i)]],
            input[[paste0("boxFLvsSSdat", i)]]))
        # potentially alter parameter
        if (selected == 1)
            ssPars[inds$ss][i] <- flPars[inds$fl][i]
    }
    # construct new flowFrame
    vals$ffControls_metsChecked <- alter_pars(ff, ssPars)
    removeModal()
})

# ------------------------------------------------------------------------------
# 3rd modal's bsButton:
# remove duplicate masses from multiplexed data
observeEvent(input$msCheckedMP, {
    inds <- duplicateMsMP()
    remove <- NULL
    for (i in seq_along(inds)) {
        choices <- mpChs()$Metal[inds[[i]]]
        selected <- input[[paste0("metalSelectionMP", i)]]
        remove <- append(remove, which(inds[[i]])[choices != selected])
    }
    ff <- fsComp()
    vals$fsComp_msChecked <- ff[, -remove]
    removeModal()
})

# ------------------------------------------------------------------------------
# 4th modal's bsButton:
# update metals of multiplexed data
observeEvent(input$metsCheckedSSvsMP, {
    inds <- matchingMsMP()
    n <- sum(inds[[1]])
    
    fs <- vals$fsComp_msChecked
    # get parameters
    mpPars <- colnames(fs)
    ssPars <- colnames(ffControls())
    for (i in seq_len(n)) {
        # determine if single-stain or multiplexed metal has been selected
        selected <- which(c(
            input[[paste0("boxSSvsMPref", i)]],
            input[[paste0("boxSSvsMPdat", i)]]))
        # potentially alter parameter
        if (selected == 1)
            mpPars[inds$mp][i] <- ssPars[inds$ss][i]
    }
    # construct new flowFrame
    vals$fsComp_metsChecked <- fsApply(fs, alter_pars, mpPars)
    removeModal()
})

# ------------------------------------------------------------------------------
# get masses and metals from measurement parameters
# of Fluidigm as reference, single-stains and multiplexed data
flChs <- reactive({
    chs <- readRDS('data/Fluidigm_default_metals.rds')
    get_ms_and_mets(chs)
})
ssChs <- reactive({
    chs <- colnames(ffControls())
    get_ms_and_mets(chs)
})
mpChs <- reactive({
    chs <- colnames(fsComp())
    get_ms_and_mets(chs)
})

# ------------------------------------------------------------------------------
# get indices of channels that match b/w Fluidigm & single-stains 
# and single-stains and multiplexed data, respectively 
matchingMsSS <- reactive({
    req(!is.null(vals$ffControls_msChecked), flChs(), ssChs())
    ms1 <- flChs()$Mass
    ms2 <- ssChs()$Mass
    ex1 <- is.na(as.numeric(ms1))
    ex2 <- is.na(as.numeric(ms2))
    ms <- intersect(ms1[!ex1], ms2[!ex2])
    inds1 <- ms1 %in% ms
    inds2 <- ms2 %in% ms
    setNames(list(inds1, inds2), c("fl", "ss"))
})
matchingMsMP <- reactive({
    req(!is.null(vals$fsComp_msChecked), ssChs(), mpChs())
    ms1 <- ssChs()$Mass
    ms2 <- mpChs()$Mass
    ex1 <- is.na(as.numeric(ms1))
    ex2 <- is.na(as.numeric(ms2))
    ms <- intersect(ms1[!ex1], ms2[!ex2])
    inds1 <- ms1 %in% ms
    inds2 <- ms2 %in% ms
    setNames(list(inds1, inds2), c("ss", "mp"))
})

# ------------------------------------------------------------------------------
# get indices of metals that differ b/w Fluidigm & single-stains 
# and single-stains and multiplexed data, respectively 
mismatchingMetsSS <- reactive({
    req(matchingMsSS())
    flMets <- flChs()$Metal[matchingMsSS()$fl]
    ssMets <- ssChs()$Metal[matchingMsSS()$ss]
    flMets != ssMets
})
mismatchingMetsMP <- reactive({
    req(matchingMsSS())
    ssMets <- ssChs()$Metal[matchingMsMP()$ss]
    mpMets <- mpChs()$Metal[matchingMsMP()$mp]
    ssMets != mpMets
})

# ------------------------------------------------------------------------------
# toggle checkboxInputs for metal selection
# ------------------------------------------------------------------------------
observeEvent(unlist(sapply(seq_len(sum(matchingMsSS()[[1]])), 
    function(i) input[[paste0("boxFLvsSSref", i)]])), {
        n <- sum(matchingMsSS()[[1]])
        sapply(seq_len(n), function(i) 
            if (input[[paste0("boxFLvsSSref", i)]])
                updateCheckboxInput(session, 
                    inputId=paste0("boxFLvsSSdat", i), 
                    value=FALSE))
    })
observeEvent(unlist(sapply(seq_len(sum(matchingMsSS()[[1]])), 
    function(i) input[[paste0("boxFLvsSSdat", i)]])), {
        n <- sum(matchingMsSS()[[1]])
        sapply(seq_len(n), function(i) 
            if (input[[paste0("boxFLvsSSdat", i)]])
                updateCheckboxInput(session, 
                    inputId=paste0("boxFLvsSSref", i), 
                    value=FALSE))
    })
observeEvent(unlist(sapply(seq_len(sum(matchingMsMP()[[1]])), 
    function(i) input[[paste0("boxSSvsMPref", i)]])), {
        n <- sum(matchingMsMP()[[1]])
        sapply(seq_len(n), function(i) 
            if (input[[paste0("boxSSvsMPref", i)]])
                updateCheckboxInput(session, 
                    inputId=paste0("boxSSvsMPdat", i), 
                    value=FALSE))
    })
observeEvent(unlist(sapply(seq_len(sum(matchingMsMP()[[1]])), 
    function(i) input[[paste0("boxSSvsMPdat", i)]])), {
        n <- sum(matchingMsMP()[[1]])
        sapply(seq_len(n), function(i) 
            if (input[[paste0("boxSSvsMPdat", i)]])
                updateCheckboxInput(session, 
                    inputId=paste0("boxSSvsMPref", i), 
                    value=FALSE))
    })

# color textInputs depending on (mis)match between metals
observeEvent(input$boxFLvsSSdat1, once=TRUE, {
    n <- length(mismatchingMetsSS())
    cols <- rep("#cce698", n)
    cols[mismatchingMetsSS()] <- "#ffcccc"
    for (i in 1:2) {
        id <- paste0("metalFLvsSS", c("ref", "dat")[i])
        sapply(seq_len(n), function(j)
            js$restyleMetal(id=paste0(id, j), col=cols[j]))
    }
})
observeEvent(input$boxSSvsMPdat1, once=TRUE, {
    n <- length(mismatchingMetsMP())
    cols <- rep("#cce698", n)
    cols[mismatchingMetsMP()] <- "#ffcccc"
    for (i in 1:2) {
        id <- paste0("metalSSvsMP", c("ref", "dat")[i])
        sapply(seq_len(n), function(j)
            js$restyleMetal(id=paste0(id, j), col=cols[j]))
    }
})