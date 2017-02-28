# load required packages
pkgs <- c("shiny", "shinydashboard", "ggplot2", "magrittr")
lapply(pkgs, require, character.only=TRUE)

# set maximum web request size to 500 MB
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
    
    source("ui/guides.R")
    source("ui/aes.R")
    source("ui/debarcoding_tab.R")
    source("ui/debarcoding_plots.R")
    source("ui/compensation_plots.R")
    source("helpers.R")
    source("scatter.R")
    
# --------------------------------------------------------------------------------------------------
    
    vals <- reactiveValues(
# DEBARCODING
        ff1 = NULL, # input flowFrame for debarcoding
        key = NULL, # debarcoding scheme
        re1 = NULL, # result assignPrelim()
        re2 = NULL, # result estCutoffs()
        re3 = NULL, # result applyCutoffs()
        ids = NULL, # current IDs with event assignments
        yp_choices  = NULL, # all barcode IDs including 0
        ep_choices  = NULL, # only IDs with events assigned
        adj_choices = NULL, # all barcode IDs excluding 0
        mhl_choices = NULL, # IDs with events assigned excluding 0
# COMPENSATION
        ff2  = NULL, # input flowFrame for compensation
        log  = NULL, # console output
        trm  = NULL, # trim value
        sm   = NULL, # spillover matrix
        ch1  = NULL, # scatter channel 1
        ch2  = NULL, # scatter channel 2
        cmp1 = NULL, # compensated flowFrame 1
        cmp2 = NULL) # compensated flowFrame 2

# ==================================================================================================
# DEBARCODING
# ==================================================================================================
    
    output$debarcoding_guide  <- renderUI({ debarcoding_guide })
    output$compensation_guide <- renderUI({ compensation_guide })
    
    output$debarcoding_sidebar_1 <- renderUI({
        if (is.null(input$fcs)) return()
        vals$ff1 <- flowCore::read.FCS(input$fcs$datapath)
        debarcoding_sidebar_1
    })

# --------------------------------------------------------------------------------------------------
# BUTTON - "Assign preliminary IDs"   
# --------------------------------------------------------------------------------------------------
    
    observeEvent(input$button_assignPrelim, {
        if (is.null(input$input_bcChs) & is.null(input$csv)) return()
        
        # get barcode key
        if (!is.null(input$input_bcChs)) {
            vals$key <- gsub("[[:alpha:][:punct:]]", "", input$input_bcChs)
        } else if (!is.null(input$csv)) {
            vals$key <- read.csv(input$csv$datapath, check.names=FALSE, row.names=1)
        }

        # assign IDs
        showNotification(h5("Assigning preliminary IDs..."), id="msg", 
                         type="message", duration=NULL, closeButton=FALSE)
        vals$re1 <- assignPrelim(x=vals$ff1, y=vals$key)
        removeNotification(id="msg")
        
        vals$adj_choices <- rownames(bc_key(vals$re1))
        vals$yp_choices  <- c(0, rownames(bc_key(vals$re1)))
        names(vals$yp_choices) <- c("All", rownames(bc_key(vals$re1)))
        vals$ids <- sort(unique(bc_ids(vals$re1)))
        vals$ep_choices  <- vals$ids
        vals$mhl_choices <- vals$ids[vals$ids != 0]
        
# --------------------------------------------------------------------------------------------------
# event plot panel
# --------------------------------------------------------------------------------------------------
        
        output$ep_panel <- renderUI ({ ep_panel(vals$ep_choices) })
        output$plot_plotEvents <- renderPlot(
            plotEvents(x = vals$re1, 
                       which = input$ep_which, 
                       n_events = as.numeric(input$n_events)))
        
        output$debarcoding_sidebar_2 <- renderUI ({ debarcoding_sidebar_2 })
        })

        observe({
            if (is.null(vals$re2)) return()
        
# --------------------------------------------------------------------------------------------------
# yield plot panel
# --------------------------------------------------------------------------------------------------
        
            output$yp_panel <- renderUI ({ yp_panel(vals$yp_choices) })
            output$plot_plotYields <- renderPlot(
                plotYields(x = vals$re2, 
                           which = input$yp_which))
            
            #tags$head(tags$style("#table_summary {white-space:nowrap}"))    
            output$table_summary <- DT::renderDataTable({ summary_tbl(vals$re2) })
        })

# --------------------------------------------------------------------------------------------------
# next / previous buttons
# --------------------------------------------------------------------------------------------------
    
    observeEvent(input$yp_prev, { 
        if (input$yp_which == 0) return()
        updateSelectInput(session, "yp_which", 
                          selected=vals$yp_choices[which(vals$yp_choices == input$yp_which)-1]) })
    
    observeEvent(input$yp_next, { 
        if (input$yp_which == vals$yp_choices[length(vals$yp_choices)]) return()
        updateSelectInput(session, "yp_which", 
                          selected=vals$yp_choices[which(vals$yp_choices == input$yp_which)+1]) })
    
    observeEvent(input$ep_prev, { 
        if (input$ep_which == vals$ep_choices[1]) return()
        updateSelectInput(session, "ep_which", 
                          selected=vals$ep_choices[which(vals$ep_choices == input$ep_which)-1]) })
    
    observeEvent(input$ep_next, { 
        if (input$ep_which == vals$ep_choices[length(vals$ep_choices)]) return()
        updateSelectInput(session, "ep_which", 
                          selected=vals$ep_choices[which(vals$ep_choices == input$ep_which)+1]) })
    
    observeEvent(input$mhl_prev, { 
        if (input$mhl_which == vals$mhl_choices[1]) return()
        updateSelectInput(session, "ep_which", 
                          selected=vals$mhl_choices[which(vals$mhl_choices == input$mhl_which)-1]) })
    
    observeEvent(input$mhl_next, { 
        if (input$mhl_which == vals$mhl_choices[length(vals$mhl_choices)]) return()
        updateSelectInput(session, "ep_which", 
                          selected=vals$mhl_choices[which(vals$mhl_choices == input$mhl_which)+1]) })
    
# --------------------------------------------------------------------------------------------------
# debarcoding tab sidebar
# --------------------------------------------------------------------------------------------------  
    
    # CHECKBOX - "Upload barcoding scheme"
    output$file_csv <- renderUI({
        if (input$box_csv)
            fileInput("csv", NULL, accept=".csv")
    })
    
    # CHECKBOX - "Select barcode channels"
    output$select_bcChs <- renderUI({
        if (input$box_bcChs)
            selectInput("input_bcChs", NULL, multiple=TRUE, size=8, selectize=FALSE,
                        choices=flowCore::colnames(vals$ff1))
    })
    
    # CHECKBOX - "Estimate separation cutoffs"
    observe({
        if (is.null(input$box_estCutoffs) ||
            input$box_estCutoffs == 0)
            return()
            vals$re2 <- estCutoffs(x = vals$re1)
    })
    
    # CHECKBOX - "Adjust population-specific cutoffs"
    # SELECT - choose bc for cutoff adjustment
    output$select_adjustCutoff <- renderUI ({
        if (input$box_indivCutoffs)
            selectInput("select_adjustCutoff", NULL,
                        choices=vals$adj_choices, selected=vals$adj_choices[1])
    })
    
    # synchronize select_adjustCutoff w/ yield plot
    observe({
        if (is.null(input$box_indivCutoffs) || 
            input$box_indivCutoffs == 0 || 
            input$yp_which == 0) 
            return() 
        updateSelectInput(session, "select_adjustCutoff", 
                          selected=input$yp_which)
    })
    
    # INPUT - enter cutoff val
    output$input_adjustCutoff <- renderUI ({
        if (input$box_indivCutoffs)
            numericInput("input_adjustCutoff", NULL,
                         min=0, max=1, step=.01,
                         value=vals$re2@sep_cutoffs[
                             input$select_adjustCutoff == rownames(vals$re2@bc_key)])
    })
    
    # BUTTON - "Adjust"
    output$button_adjustCutoff <- renderUI ({
        if (input$box_indivCutoffs) 
            actionButton("button_adjustCutoff", "Adjust", style=ylw_button)
    })
    
    observeEvent(input$button_adjustCutoff, {
        ind <- rownames(vals$re2@bc_key) == input$select_adjustCutoff
        vals$re2@sep_cutoffs[ind] <- input$input_adjustCutoff
    })
    
    
    # CHECKBOX - "Enter global separation cutoff"
    # input_globalCutoff
    output$input_globalCutoff <- renderUI ({
        if (input$box_globalCutoff)
            numericInput("input_globalCutoff", NULL, width="100%",
                         value=.2, min=0, max=1, step=.01)
    })
    
    # button_globalCutoff
    output$button_globalCutoff <- renderUI ({
        if (input$box_globalCutoff)
            actionButton("button_globalCutoff", "Apply", style=ylw_button)
    })
    
    observeEvent(input$button_globalCutoff, {
        sep_cutoffs(vals$re2) <- input$input_globalCutoff
    })
    
# --------------------------------------------------------------------------------------------------
# synchronize plots
# -------------------------------------------------------------------------------------------------- 
    
    observe({
        if (is.null(input$yp_which) || 
            input$yp_which == 0) return() 
        updateSelectInput(session, "ep_which",  selected=input$yp_which)
        updateSelectInput(session, "mhl_which", selected=input$yp_which)
    })
    
    observe({
        if (is.null(input$ep_which) ||
            input$ep_which == 0) return() 
        updateSelectInput(session, "yp_which",  selected=input$ep_which)
        updateSelectInput(session, "mhl_which", selected=input$ep_which)
    })
    
    observe({
        if (is.null(input$mhl_which)) return()
        updateSelectInput(session, "yp_which",  selected=input$mhl_which)
        updateSelectInput(session, "ep_which",  selected=input$mhl_which)
    })
    
# --------------------------------------------------------------------------------------------------
# toggle checkboxes
# -------------------------------------------------------------------------------------------------- 
    
    observe({
        if (is.null(input$box_csv) || 
            input$box_csv == 0) 
            return()
        updateCheckboxInput(session, "box_bcChs", value=FALSE)
    })
    
    observe({
        if (is.null(input$box_bcChs) || 
            input$box_bcChs == 0) 
            return()
        updateCheckboxInput(session, "box_csv", value=FALSE)
    })
    
    observe({
        if (is.null(input$box_estCutoffs) || 
            input$box_estCutoffs == 0) 
            return()
        updateCheckboxInput(session, "box_indivCutoffs", value=FALSE)
        updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
    })
    
    observe({
        if (is.null(input$box_indivCutoffs) || 
            input$box_indivCutoffs == 0) 
            return()
        updateCheckboxInput(session, "box_estCutoffs",   value=FALSE)
        updateCheckboxInput(session, "box_globalCutoff", value=FALSE)
    })
    
    observe({
        if (is.null(input$box_globalCutoff) || 
            input$box_globalCutoff == 0) 
            return()
        updateCheckboxInput(session, "box_estCutoffs",   value=FALSE)
        updateCheckboxInput(session, "box_indivCutoffs", value=FALSE)
    })

# --------------------------------------------------------------------------------------------------
# BUTTON - "Apply cutoffs"
# --------------------------------------------------------------------------------------------------
    
    observeEvent(input$button_applyCutoffs, {
        
        showNotification(h6("Applying thresholds..."), id="msg", 
                         type="message", duration=NULL, closeButton=FALSE)
        vals$re3 <- applyCutoffs(x = vals$re2,
                                 mhl_cutoff = input$input_mhlCutoff)
        removeNotification(id="msg")
        
        # update choices
        vals$ids <- sort(unique(bc_ids(vals$re3)))
        vals$ep_choices  <- vals$ids
        vals$mhl_choices <- vals$ids[vals$ids != 0]
        
        # update event panel
        output$plot_plotEvents <- renderPlot(
            plotEvents(x = vals$re3,
                       which = input$ep_which,
                       n_events = as.numeric(input$n_events)))
        
        # mahal plot panel
        output$mhl_panel <- renderUI ({ mhl_panel(vals$mhl_choices) })
        output$plot_plotMahal <- renderPlot( 
            plotMahal(x = vals$re3, 
                      which = input$mhl_which,
                      cofactor = input$input_mhlCofactor))
    })

# --------------------------------------------------------------------------------------------------
# checkboxes "Use medians" | "Use trim estimate" | "Enter trim value"
# --------------------------------------------------------------------------------------------------
    
    # toggle checkboxes
    observe({
        if (is.null(input$box_useMedians) || 
            input$box_useMedians == 0) 
            return()
        updateCheckboxInput(session, "box_estTrim",   value=FALSE)
        updateCheckboxInput(session, "box_enterTrim", value=FALSE)
    })
    
    observe({
        if (is.null(input$box_estTrim) || 
            input$box_estTrim == 0) 
            return()
        updateCheckboxInput(session, "box_useMedians", value=FALSE)
        updateCheckboxInput(session, "box_enterTrim",  value=FALSE)
    })
    
    observe({
        if (is.null(input$box_enterTrim) || 
            input$box_enterTrim == 0) 
            return()
        updateCheckboxInput(session, "box_useMedians", value=FALSE)
        updateCheckboxInput(session, "box_estTrim",    value=FALSE)
    })
    
    # checkbox "Use medians"
    observe({
        if (is.null(input$box_useMedians) || input$box_useMedians == 0) return()
        showNotification(h5("Estimating spillover..."), id="msg", 
                         type="message", duration=NULL, closeButton=FALSE)
        vals$sm  <- computeSpillmat(x=vals$re3, method="median")
        removeNotification(id="msg")
    })

    # checkboxe "Use trim estimate"
    observe({
        if (is.null(input$box_estTrim) || input$box_estTrim == 0) return()
        showNotification(h5("Estimating optimal trim value..."), id="msg", 
                         type="message", duration=NULL, closeButton=FALSE)
        vals$trm <- estTrim(x=vals$re3)
        removeNotification(id="msg")
        showNotification(h5("Estimating spillover..."), id="msg", 
                         type="message", duration=NULL, closeButton=FALSE)
        vals$sm  <- computeSpillmat(x=vals$re3, trim=vals$trm)
        removeNotification(id="msg")
    })
    
    # checkbox "Enter trim value"
    output$text_enterTrim <- renderUI({
        if (input$box_enterTrim) 
            helpText("Note that a trim value of 0.5 is equal to using medians.")
    })
    output$input_enterTrim <- renderUI ({
        if (input$box_enterTrim)
        numericInput("input_enterTrim", label=NULL,
                     value=NULL, min=.01, max=.5, step=.01) 
    })
    output$button_enterTrim <- renderUI({
        if (input$box_enterTrim) 
            actionButton("button_enterTrim", "Confirm", style=ylw_button)
    })
    observeEvent(input$button_enterTrim, {
        if (is.null(input$input_enterTrim)) return()
        showNotification(h5("Estimating spillover..."), id="msg", 
                         type="message", duration=NULL, closeButton=FALSE)
        vals$trm <- input$input_enterTrim
        vals$sm  <- computeSpillmat(x=vals$re3, trim=vals$trm)
        removeNotification(id="msg")
    })
    
####################################################################################################
# COMPENSATION
####################################################################################################
    
# --------------------------------------------------------------------------------------------------
# sidebar
    
    observe({
        if (is.null(input$fcs2) ||
            is.null(vals$re3)) 
            return()
        vals$ff2 <- flowCore::read.FCS(input$fcs2$datapath)    
        output$compensation_sidebar_1 <- renderUI({ compensation_sidebar_1 })
    })
    
    output$input_upldSM <- renderUI({
        if (input$box_upldSM)
            fileInput("input_SM", NULL, accept=".csv")
    })
    
    observe({
        if (!is.null(input$input_SM))
            vals$sm <- read.csv(input$input_SM$datapath, check.names=FALSE, row.names=1)
    })
    
    observe({
        if (is.null(input$box_estSM) || 
            input$box_estSM == 0 || 
            is.null(vals$re3)) 
            return()
        showNotification(h5("Estimating spillover..."), id="msg", 
                         type="message", duration=NULL, closeButton=FALSE)
        vals$sm <- computeSpillmat(x = vals$re3)
        removeNotification(id="msg")
    })
    
    output$text_compCytof <- renderUI({
        if (is.null(input$box_estSM) ||
            input$box_estSM == 0)
            return()
        tagList(textOutput("text_compCytof_1"), 
                verbatimTextOutput("text_compCytof_2"))
    })

# --------------------------------------------------------------------------------------------------
# toggle checkboxes
    
    observe({
        if (is.null(input$box_upldSM) ||
            input$box_upldSM == 0)
            return()
        updateCheckboxInput(session, "box_estSM", value=FALSE)
    })
    
    observe({
        if (is.null(input$box_estSM) ||
            input$box_estSM == 0)
            return()
        updateCheckboxInput(session, "box_upldSM", value=FALSE)
    })
    
# -------------------------------------------------------------------------------------------------- 
    
    observe({
        if (is.null(vals$ff1) || 
            is.null(input$input_bcChs))
            return()
        output$panel_estTrim      <- renderUI({ panel_estTrim      })
        output$panel_plotSpillmat <- renderUI({ panel_plotSpillmat })
        output$panel_plotScatter  <- renderUI({ panel_plotScatter  }) 
        output$panel_scatters     <- renderUI({ 
            panel_scatters(flowCore::colnames(vals$ff1),
                input$input_bcChs[1], input$input_bcChs[2])
        })
    })
    
    observe({
        if (is.null(vals$sm)) return()
        output$text_compCytof_1 <- renderText(
            "WARNING: Compensation is likely to be inaccurate. 
            Spill values for the following interactions have not been estimated:")
        output$text_compCytof_2 <- renderPrint(
            cat(capture.output(vals$cmp1 <- compCytof(x=vals$ff1, y=vals$sm)), sep="\n"))
        vals$cmp2 <- compCytof(x=vals$ff2, y=vals$sm)
        
    })

    output$plot_plotSpillmat <- renderPlot({
        if (is.null(vals$sm)) return()
        input$button_newSpill
        plotSpillmat(bc_ms=vals$key, SM=isolate(vals$sm))
    })
    
    output$plot_plotScatter <- renderPlot({
        if (is.null(vals$sm)) return()
        input$button_newSpill
        plotScatter(x=vals$re3, SM=isolate(vals$sm))
    })
    
    output$compensation_sidebar_2 <- renderUI({
        if (is.null(vals$cmp1) || is.null(vals$cmp2) || is.null(vals$sm)) return()
        compensation_sidebar_2
    })

# --------------------------------------------------------------------------------------------------
# scatters
# --------------------------------------------------------------------------------------------------
    
    output$text_spill <- renderText({
        input$button_view
        paste0(sprintf("%.3f", 100*vals$sm[isolate(input$input_scatterCh1), 
                                           isolate(input$input_scatterCh2)]), "%")
    })
    
    # adjust spill
    observeEvent(input$button_newSpill, {
        if (is.null(input$input_newSpill)) 
            return()
        vals$sm[input$input_scatterCh1, 
                input$input_scatterCh2] <- input$input_newSpill/100
    })
    
    # top-left scatter (ff1 uncompensated)
    output$plot_scatter1 <- renderPlot({ 
        input$button_view
        scatter(ff=vals$ff1, 
                which=c(isolate(input$input_scatterCh1), 
                        isolate(input$input_scatterCh2)),
                cofactor=isolate(input$input_cofactor))
    })
    output$text_info1 <- renderText({ 
        text_info(vals$ff1, isolate(input$input_cofactor), input$rect1, 
            isolate(input$input_scatterCh1), isolate(input$input_scatterCh2))
    })   
    
    # top-right scatter (ff1 compensated)
    output$plot_scatter2 <- renderPlot({
        if (is.null(vals$cmp1)) return()
        (input$button_view || input$button_newSpill)
        scatter(ff=isolate(vals$cmp1),
                which=c(isolate(input$input_scatterCh1),
                        isolate(input$input_scatterCh2)),
                cofactor=isolate(input$input_cofactor))
    })
    output$text_info2 <- renderText({ 
        input$button_view 
        text_info(isolate(vals$cmp1), isolate(input$input_cofactor), input$rect2, 
            isolate(input$input_scatterCh1), isolate(input$input_scatterCh2))
    })
    
    # bottom-left scatter (ff2 uncompensated)
    output$plot_scatter3 <- renderPlot({ 
        input$button_view 
        scatter(ff=vals$ff2,
                which=c(isolate(input$input_scatterCh1), 
                        isolate(input$input_scatterCh2)),
                cofactor=isolate(input$input_cofactor))
    })
    output$text_info3 <- renderText({ 
        text_info(vals$ff2, isolate(input$input_cofactor), input$rect3,
                  isolate(input$input_scatterCh1), isolate(input$input_scatterCh2))
    }) 
    
    # bottom-left scatter (ff2 uncompensated)
    output$plot_scatter4 <- renderPlot({ 
        if (is.null(vals$cmp2)) return()
        (input$button_view || input$button_newSpill)
        scatter(ff=isolate(vals$cmp2),
                which=c(isolate(input$input_scatterCh1), 
                        isolate(input$input_scatterCh2)),
                cofactor=isolate(input$input_cofactor))
    })
    output$text_info4 <- renderText({ 
        input$button_view 
        text_info(isolate(vals$cmp2), isolate(input$input_cofactor), input$rect4, 
            isolate(input$input_scatterCh1), isolate(input$input_scatterCh2))
    })
    
# --------------------------------------------------------------------------------------------------
# download buttons
# --------------------------------------------------------------------------------------------------
    
    output$dwnld_fcs  <- downloadHandler(filename = function()     { "fcs.zip" },
                                         content  = function(file) { tmpdir <- tempdir(); setwd(tmpdir)
                                         outFCS(x = vals$re3, out_path = tmpdir) 
                                         zip(zipfile=file, files=c(paste0(ep, ".fcs"))) }, 
                                         contentType = "application/zip")
    
    output$dwnld_yep   <- downloadHandler(filename = function()     { "yield_event_plots.zip" },
                                          content  = function(file) { tmpdir <- tempdir(); setwd(tmpdir)
                                          plotYields(x = vals$re3, which = vals$yp_choices, out_path = tmpdir)
                                          plotEvents(x = vals$re3, which = vals$ep_choices, out_path = tmpdir, n_events = as.numeric(input$n_events))
                                          zip(zipfile=file, contentType = "application/zip", files=paste0(c("summary_yield_plot", "yield_plot", "event_plot"), ".pdf")) })
    
    output$dwnld_comped_1 <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs), "-comped.fcs") },
                                             content  = function(file) { flowCore::write.FCS(vals$cmp1, file) })
    output$dwnld_comped_2 <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs2),"-comped.fcs") },
                                             content  = function(file) { flowCore::write.FCS(vals$cmp2, file) })
    output$dwnld_spillMat <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs), "-spillMat.csv") },
                                             content  = function(file) { write.csv(vals$sm, file) })
    
})



