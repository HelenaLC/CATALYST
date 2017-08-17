# ==============================================================================
# CONCATENATION
# ==============================================================================

# read input FCS
ffsConcat <- reactive({
    req(input$fcsConcat)
    n <- nrow(input$fcsConcat)
    # check validity of input FCS files
    valid <- check_FCS_fileInput(input$fcsConcat, n)
    if (!valid) return()
    ffs <- lapply(seq_len(n), function(i)
        flowCore::read.FCS(
            filename=input$fcsConcat[[i, "datapath"]],
            transformation=FALSE,
            truncate_max_range=FALSE))
    # check validity of input flowFrames
    nPars <- sapply(ffs, function(ff) ncol(ff))
    invalid <- length(unique(nPars)) != 1
    if (invalid) {
        showNotification(
            "Input flowFrames do not contain 
            the same number of measurement parameters.",
            duration=NULL, type="error")
        return()
    }
    pars <- sapply(ffs, function(ff) colnames(ff))
    unique <- apply(pars, 1, function(par) length(unique(par)) == 1)
    invalid <- sum(unique) != nPars[1]
    if (invalid) {
        showNotification(
            "Input flowFrames do not contain identical stains.",
            duration=NULL, type="error")
        return()
    }
    as(ffs, "flowSet")
})

# render shinydashboard::box "File editing"
output$editFcsConcat <- renderUI(editFcsConcat(ffsConcat()[[1]]))

# render checkboxInputs "Add file number channel" & "Order by time"
# & downloadButton "Merge files"
output$concatOutput <- renderUI({
    req(ffsConcat())
    nm <- input$fcsConcat$name[1]
    outNm <- gsub("(_)*([[:digit:]]*)(.fcs)$", 
        "_concat.fcs", nm, ignore.case=TRUE)
    concatOutput(outNm)
})

# bsButton "Go to normalization": propagate data 
# & disable FCS fileInput from normalization tab
observeEvent(input$goToNorm, {
    vals$keepDataConcat <- TRUE
    shinyjs::disable(id="fcsNorm")
    updateTabItems(session, inputId="tabs", selected="normalization")
})

# ------------------------------------------------------------------------------
# download handler
# ------------------------------------------------------------------------------

ffConcat <- reactive({
    req(isTRUE(vals$keepDataConcat))
    nPars <- ncol(ffsConcat()[[1]])
    pars <- sapply(seq_len(nPars), function(i) input[[paste0("par", i)]])
    desc <- sapply(seq_len(nPars), function(i) input[[paste0("des", i)]])
    CATALYST::concatFCS(ffsConcat(), 
        by_time=input$orderByTime,
        file_num=input$addFileNo, 
        pars=pars, 
        desc=desc)
})

# get output file name
smplNmConcat <- reactive({
    nm <- input$concatOutputFileNm
    req(nm)
    # assure suffix is ".fcs"
    nm <- gsub("([.][[:alpha:]]*$)", "", nm, ignore.case=TRUE)
    if (nchar(nm) == 0)
        nm <- paste0(format(Sys.Date(), "%y%m%d"), "_concat")
    paste0(nm, ".fcs")
})

output$dwnldConcat <- downloadHandler(
    filename=function() {
        smplNmConcat()
    }, 
    content=function(file) { 
        fs <- as(ffsConcat(), "flowSet")
        n <- length(fs)
        # order by time
        if (input$orderByTime) {
            bts <- keyword(fs, "$BTIM")
            o <- order(bts)
            fs <- fs[o]
        }
        nPars <- ncol(fs[[1]])
        nEvents <- as.numeric(keyword(fs, "$TOT"))
        # concatenate
        es <- fsApply(fs, exprs)
        timeCol <- grep("time", colnames(fs), ignore.case=TRUE)
        start <- c(1, cumsum(nEvents)+1)
        end <- start[-1]-1
        for (i in seq_along(fs)[-1]) {
            inds <- start[i]:end[i]
            es[inds, timeCol] <- es[inds, timeCol]+es[end[i-1], timeCol]
        }
        # ······································································
        # construct new flowFrame
        # (mandatory keywords are set by flowCore::write.FCS)
        # ······································································
        # get parameters and descriptions from textInputs
        pars <- sapply(seq_len(nPars), function(i) input[[paste0("par", i)]])
        desc <- sapply(seq_len(nPars), function(i) input[[paste0("des", i)]])
        # get descriptions and modify some keywords
        d <- description(fs[[1]])
        d <- d[!names(d) %in% c("$FIL", "FILENAME", "ORIGINALGUID")]
        d$`GUID` <- gsub("_\\d+.fcs", "_concat.fcs", 
            description(fs[[1]])$ORIGINALGUID, ignore.case=TRUE)
        d$`FILE LIST` <- paste(paste0(seq_len(n), 
            "-", keyword(fs, "ORIGINALGUID")), collapse=",")
        d$`$BTIM` <- description(fs[[1]])$`$BTIM`
        d$`$ETIM` <- description(fs[[n]])$`$ETIM`
        d$`COM` <- "FCS files concatenated by CATALYST"
        d[paste0("$P", seq_len(nPars), "N")] <- pars
        d[paste0("$P", seq_len(nPars), "S")] <- desc
        # get ranges
        rngs <- parameters(fs[[1]])$range
        mins <- parameters(fs[[1]])$minRange
        maxs <- parameters(fs[[1]])$maxRange
        # add file number channel
        if (input$addFileNo) {
            nPars <- nPars + 1
            rngs <- append(rngs, -1)
            mins <- append(mins, 0)
            maxs <- append(maxs, 1)
            pars <- append(pars, "FileNum")
            desc <- append(desc, "File Number")
            d[paste0("$P", nPars, "B")] <- 32
            d[paste0("$P", nPars, "E")] <- 0
            d[paste0("$P", nPars, "N")] <- "FileNum"
            d[paste0("$P", nPars, "R")] <- 0
            d[paste0("$P", nPars, "S")] <- "File Number"
            es <- cbind(es, rep(seq_len(n)[o], keyword(fs, "$TOT")))
        }
        # create AnnotatedDataFrame of parameters
        df <- data.frame(list(
            name=pars, 
            desc=desc,
            range=rngs,
            minRange=mins,
            maxRange=maxs))
        p <- Biobase::AnnotatedDataFrame(df)
        rownames(p) <- paste0("$P", seq_len(nPars))
        # add file number channel
        # make flowFrame
        colnames(es) <- pars
        ff <- flowFrame(exprs=es, parameters=p, description=d)
        identifier(ff) <- smplNmConcat()
        #write.FCS(ff, file)
        suppressWarnings(write.FCS(ff, file))
    }
)