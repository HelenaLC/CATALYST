# load required packages
pkgs <- c("shiny")
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
                         type="warning", duration=NULL, closeButton=FALSE)
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
                         type="warning", duration=NULL, closeButton=FALSE)
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
# use mean / estimate/enter trim value
# --------------------------------------------------------------------------------------------------
    
    observe({
        if (is.null(input$use_mean) || 
            input$use_mean == 0) 
            return()
        updateCheckboxInput(session, "est_trim",   value=FALSE)
        updateCheckboxInput(session, "enter_trim", value=FALSE)
    })
    
    observe({
        if (is.null(input$est_trim) || 
            input$est_trim == 0) 
            return()
        updateCheckboxInput(session, "use_mean",   value=FALSE)
        updateCheckboxInput(session, "enter_trim", value=FALSE)
    })
    
    observe({
        if (is.null(input$enter_trim) || 
            input$enter_trim == 0) 
            return()
        updateCheckboxInput(session, "use_mean", value=FALSE)
        updateCheckboxInput(session, "est_trim", value=FALSE)
    })
    
    output$enter_trim <- renderUI ({
        if (input$enter_trim) 
            numericInput("input_trim", label=NULL, width="32%",
                         value=.05, min=.01, max=.5, step=.01) })
    
####################################################################################################
# COMPENSATION
####################################################################################################
    
    observe({
        if (is.null(input$fcs2) ||
            is.null(vals$re3)) 
            return()
        
        # read input flowFrame
        vals$ff2 <- flowCore::read.FCS(input$fcs2$datapath)
        
        # render sidebar
        output$compensation_sidebar_1 <- renderUI({ compensation_sidebar_1 })

        # render plot tabs
        output$panel_estTrim      <- renderUI({ panel_estTrim      })
        output$panel_plotSpillmat <- renderUI({ panel_plotSpillmat })
        output$panel_plotScatter  <- renderUI({ panel_plotScatter  }) 
        output$panel_scatters     <- renderUI({ 
            panel_scatters(flowCore::colnames(vals$ff1),
                input$input_bcChs[1], input$input_bcChs[2])
        })
        
        # estTrim()
        output$plot_estTrim <- renderPlot(
            (estTrim(x = vals$re3)))
        
        # computeSpillmat()
        vals$sm <- computeSpillmat(x = vals$re3)
        
        # plotScatter()
        output$plot_plotScatter <- renderPlot(
            plotScatter(x = vals$re3,
                        SM = vals$sm)
        )

        # spillover matrix heat map
        output$plot_plotSpillmat <- renderPlot(
            plotSpillmat(bc_ms = vals$key,
                         SM = vals$sm)
        )
        
        # compensate
        output$text_compCytof_1 <- renderText({
                "WARNING: Compensation is likely to be inaccurate. 
                 Spill values for the following interactions have not been estimated:"
        })
        output$text_compCytof_2 <- renderPrint({
            cat(capture.output(vals$cmp1 <- compCytof(x = vals$ff1, y = vals$sm)), sep="\n")
        })
        vals$cmp2 <- compCytof(x = vals$ff2, y = vals$sm)
    })

# --------------------------------------------------------------------------------------------------
# scatters
# --------------------------------------------------------------------------------------------------
    
    output$text_spill <- renderText({
        input$button_view
        paste0(sprintf("%.3f", 100*vals$sm[isolate(input$input_scatterCh1), 
                                    isolate(input$input_scatterCh2)]), "%")
    })
    
    # top-left scatter (ff1 uncompensated)
    output$plot_scatter1 <- renderPlot({
        input$button_view
        scatter(ff = vals$ff1,
                which = c(isolate(input$input_scatterCh1),
                          isolate(input$input_scatterCh2)),
                cofactor = isolate(input$input_cofactor))
    })
    output$text_info1 <- renderText({ 
        text_info(vals$ff1, input$input_cofactor, input$rect1, 
            isolate(input$input_scatterCh1), isolate(input$input_scatterCh2))
    })   
    
    # top-right scatter (ff1 compensated)
    output$plot_scatter2 <- renderPlot({
        input$button_view
        scatter(ff = vals$cmp1,
                which = c(isolate(input$input_scatterCh1),
                          isolate(input$input_scatterCh2)),
                cofactor = isolate(input$input_cofactor))
    })
    output$text_info2 <- renderText({ 
        input$button_view 
        text_info(vals$cmp1, input$input_cofactor, input$rect2, 
            isolate(input$input_scatterCh1), isolate(input$input_scatterCh2))
    })
    
    # bottom-left scatter (ff2 uncompensated)
    output$plot_scatter3 <- renderPlot({ 
        input$button_view 
        scatter(ff = vals$ff2,
                which = c(isolate(input$input_scatterCh1), 
                          isolate(input$input_scatterCh2)),
                cofactor = isolate(input$input_cofactor))
    })
    output$text_info3 <- renderText({ 
        text_info(vals$ff2, input$input_cofactor, input$rect3,
                  isolate(input$input_scatterCh1), isolate(input$input_scatterCh2))
    }) 
    
    # bottom-left scatter (ff2 uncompensated)
    output$plot_scatter4 <- renderPlot({ 
        input$button_view 
        scatter(ff = vals$cmp2,
                which = c(isolate(input$input_scatterCh1), 
                          isolate(input$input_scatterCh2)),
                cofactor = isolate(input$input_cofactor))
    })
    output$text_info4 <- renderText({ 
        input$button_view 
        text_info(vals$cmp2, input$input_cofactor, input$rect4, 
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
    
    output$dwnld_1  <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs), "-comped.fcs") },
                                       content  = function(file) { flowCore::write.FCS(ffc,  file) })
    output$dwnld_2  <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs2),"-comped.fcs") },
                                       content  = function(file) { flowCore::write.FCS(ff2c, file) })
    output$dwnld_cm <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs), "-sm.csv") },
                                       content  = function(file) { write.csv(CM, file) })
    
})



