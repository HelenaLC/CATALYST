# load required packages
library(magrittr)

# set maximum web request size to 500 MB
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
    
    source("helpers-normalization_tmp.R")
    library(shinyjs)
    library(plotly)
    library(grid)
    library(gridExtra)
    library(CATALYST)
    library(ggplot2)
    source("helpers.R")

# guides -----   
    source("ui/guides.R")
    output$debarcoding_guide   <- renderUI(debarcoding_guide)
    output$compensation_guide  <- renderUI(compensation_guide)
    
# server files -----
    source("server/normalization-server.R", local=TRUE)
    source("server/debarcoding-server.R", local=TRUE)
    source("server/compensation-server.R", local=TRUE)
# plot files -----    
    source("ui/normalization_plots.R")
    source("ui/debarcoding_plots.R")
    source("ui/compensation_plots.R")

# ------------------------------------------------------------------------------
    
    vals <- reactiveValues(
# ----- NORMALIZATION -----
        ff_norm = NULL,     # input FCS for normalization
        ffsNorm = list(),   # input flowFrame(s) to be normalized
        ffsNormTo = NULL, # flowFrame(s) to use as baseline
        ffsNormed = list(NULL), # normalized flowFrame(s)
        beadInds = list(NULL),  # logical vector for each flowFrame
        beads = NULL,       # beads for normalization
        beadCols = NULL,    # bead columns
        dnaCols = NULL,     # DNA columns
        timeCol = NULL,     # time column
        smpl = 1,           # index of currently selected sample
# ----- DEBARCODING -----
        ff1 = NULL, # input flowFrame for debarcoding
        key = NULL, # debarcoding scheme
        re0 = NULL, # result assignPrelim()
        re1 = NULL, # result estCutoffs()
        re2 = NULL, # result applyCutoffs()
        mhl = 30,   # mhl_cutoff
        ids = NULL, # current IDs with event assignments
        yp_choices  = NULL, # all barcode IDs including 0
        ep_choices  = NULL, # only IDs with events assigned
        adj_choices = NULL, # all barcode IDs excluding 0
        mhl_choices = NULL, # IDs with events assigned excluding 0
        nms = NULL,  # names to use for output FCS files
# ----- COMPENSATION -----
        ff2  = NULL, # input flowFrame for compensation
        log  = NULL, # console output
        trm  = NULL, # trim value
        sm   = NULL, # spillover matrix
        ch1  = NULL, # scatter channel 1
        ch2  = NULL, # scatter channel 2
        cmp1 = NULL, # compensated flowFrame 1
        cmp2 = NULL) # compensated flowFrame 2

# ==============================================================================
# NORMALIZATION
# ==============================================================================
    
    output$select_customBeads <- renderUI({
        if (input$select_beads == "custom") 
            fluidRow(column(12,
                div(style="display:inline-block; width:80%", 
                    selectizeInput("input_customBeads", NULL, multiple=TRUE, 
                        choices=setNames(as.list(flowCore::colnames(vals$ffNormTo)), 
                            paste0(flowCore::parameters(vals$ffNormTo)$desc, 
                                " [", flowCore::colnames(vals$ffNormTo), "]")))),
                div(style="display:inline-block; vertical-align:top; width:19%", 
                    actionButton("button_customBeads", "Gate", icon("object-ungroup"), width="100%"))
            ))
    })
    
    observeEvent(input$button_customBeads, {
        if (is.null(input$input_customBeads)) return()
        chs <- flowCore::colnames(vals$ffNormTo)
        ms <- gsub("[[:alpha:][:punct:]]", "", chs)
        beadMs <- gsub("[[:alpha:][:punct:]]", "", unlist(input$input_customBeads))
        if (sum(beadMs %in% ms) != length(beadMs)) {
            showNotification("Not all bead channels are valid.", 
                             type="error", closeButton=FALSE)
            return()
        }
        vals$beads <- as.numeric(beadMs)
    })
    observe({
        x <- input$select_beads
        if (is.null(x)) return()
        if (x == "dvs" || x == "beta")
            vals$beads <- input$select_beads
    })
     observe({
        if (is.null(vals$beads)) return()
        chs <- flowCore::colnames(c(vals$ffsNorm, vals$ffsNormTo)[[vals$smpl]])
        vals$beadCols <- get_bead_cols(chs, vals$beads)
        vals$dnaCols  <- grep("Ir191|Ir193", chs, ignore.case=TRUE)
        vals$timeCol <- grep("time", chs, ignore.case=TRUE)
        output$box4 <- renderUI(box4)
        output$box_beadGating <- renderUI(
            box_beadGating(c(input$fcsNorm$name, input$fcsNormTo$name)))
        js$collapse("box_3")
    })
    
# ==============================================================================
# DEBARCODING
# ==============================================================================
    
    output$debarcoding_sidebar_1 <- renderUI({
        if (is.null(input$fcs)) return()
        vals$ff1 <- flowCore::read.FCS(input$fcs$datapath)
        debarcoding_sidebar_1
    })

# ------------------------------------------------------------------------------
# BUTTON: "Assign preliminary IDs"   
# ------------------------------------------------------------------------------
    
    observeEvent(input$button_assignPrelim, {
        if (is.null(input$csv) & is.null(input$input_bcChs)) {
            showNotification("Please upload a barcoding scheme or select 
                single-positive channels.", type="error", closeButton=FALSE)
            return()
        }
        
        # get barcode key
        if (!is.null(input$input_bcChs)) {
            vals$key <- gsub("[[:alpha:][:punct:]]", "", input$input_bcChs)
        } else if (!is.null(input$csv)) {
            vals$key <- read.csv(input$csv$datapath, 
                check.names=FALSE, row.names=1)
        }

        # assign IDs
        showNotification(h5("Assigning preliminary IDs..."), id="msg", 
            type="message", duration=NULL, closeButton=FALSE)
        vals$re0 <- CATALYST::assignPrelim(x=vals$ff1, y=vals$key)
        vals$re1 <- CATALYST::estCutoffs(x=vals$re0)
        removeNotification(id="msg")
        
        # inputSelect choices for yield plot and cutoff adjustment
        vals$adj_choices <- rownames(bc_key(vals$re0))
        vals$yp_choices  <- c(0, rownames(bc_key(vals$re0)))
        names(vals$yp_choices) <- c("All", rownames(bc_key(vals$re0)))

        # render yield, event and mahal plot panel
        output$yp_panel <- renderUI(yp_panel(vals$yp_choices))
        output$ep_panel <- renderUI(ep_panel(vals$ep_choices))
        output$mhl_panel <- renderUI(mhl_panel(vals$mhl_choices))
        
        output$debarcoding_sidebar_2 <- renderUI ({ debarcoding_sidebar_2 })
    })
    
    observeEvent(input$button_mhlCutoff, { 
        vals$mhl <- input$slider_mhlCutoff 
    })
    
    # apply cutoffs if deconvolution parameters change
    observe({
        if (is.null(vals$re1)) return()
        vals$re2 <- CATALYST::applyCutoffs(vals$re1, vals$mhl)
    })
    
    # estimate trim if dbFrame or sequence change
    observe({
        if (is.null(input$input_bcChs) || is.null(vals$re2)) return()
        output$plot_estTrim <- renderPlotly({
            CATALYST::estTrim(
                x=vals$re2,
                min=input$estTrim_min,
                max=input$estTrim_max,
                step=input$estTrim_step
            )
        })
    })
    
# ==============================================================================
# COMPENSATION
# ==============================================================================
    
# ------------------------------------------------------------------------------
# checkboxes "Enter trim value" & "Use medians"
# ------------------------------------------------------------------------------
    
    observeEvent(input$button_estTrim, {
        if (is.null(vals$re2)) return()
        l <- length(seq(
            min <- input$estTrim_min, 
            max <- input$estTrim_max, 
            step <- input$estTrim_step))
        showNotification(h5(paste("Estimating spill for", l, "trim values...")), id="msg", 
                         type="message", duration=NULL, closeButton=FALSE)
        output$plot_estTrim <- renderPlotly({ 
            CATALYST::estTrim(
                x=vals$re2, 
                min=min, 
                max=max, 
                step=step)
        })
        removeNotification(id="msg")
    })
    
    observeEvent(input$button_enterTrim, {
        if (is.null(input$input_enterTrim)) return()
        showNotification(h5("Estimating spillover..."), id="msg", 
                         type="message", duration=NULL, closeButton=FALSE)
        vals$trm <- input$input_enterTrim
        vals$sm  <- CATALYST::computeSpillmat(x=vals$re2, trim=vals$trm)
        removeNotification(id="msg")
    })
    
    observe({
        if (is.null(input$fcs2))
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
            vals$sm <- as.matrix(read.csv(input$input_SM$datapath, check.names=FALSE, row.names=1))
    })
    
    observe({
        if (is.null(input$box_estSM) || 
            input$box_estSM == 0 || 
            is.null(vals$re2)) 
            return()
        showNotification(h5("Estimating spillover..."), id="msg", 
            type="message", duration=NULL, closeButton=FALSE)
        vals$sm <- CATALYST::computeSpillmat(x=vals$re2)
        removeNotification(id="msg")
    })
    
    output$text_compCytof <- renderUI({
        if (is.null(input$box_estSM) ||
            input$box_estSM == 0)
            return()
        tagList(
            textOutput("text_compCytof_1"), 
            verbatimTextOutput("text_compCytof_2"))
    })

    output$upldNms <- renderUI({
        if (input$box_upldNms)
        tagList(
            helpText("A table with 2 columns: 
                Sample IDs and the desired file names."),
            fileInput("input_upldNms", NULL, accept=".csv"))
    })

    observe({
        if (is.null(input$input_upldNms) || input$box_upldNms == 0) return()
        vals$nms <- read.csv(input$input_upldNms$datapath, header=FALSE)
    })

# ------------------------------------------------------------------------------
# toggle checkboxes
    
    observe({
        x <- input$box_upldSM
        if (is.null(x) || x == 0) return()
        updateCheckboxInput(session, "box_estSM", value=FALSE)
    })
    
    observe({
        x <- input$box_estSM
        if (is.null(x) || x == 0) return()
        updateCheckboxInput(session, "box_upldSM", value=FALSE)
    })
    
    observe({
        x <- input$box_IDsAsNms
        if (is.null(x) || x == 0) return()
        updateCheckboxInput(session, "box_upldNms", value=FALSE)
    })
    
    observe({
        x <- input$box_upldNms
        if (is.null(x) || x == 0) return()
        updateCheckboxInput(session, "box_IDsAsNms", value=FALSE)
    })
    
# ------------------------------------------------------------------------------ 
    
    observe({
        if (is.null(vals$ff1) || 
            is.null(input$input_bcChs))
            return()
        output$panel_estTrim      <- renderUI({ panel_estTrim      })
        output$panel_plotSpillmat <- renderUI({ panel_plotSpillmat })
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
            cat(capture.output(
                vals$cmp1 <- compCytof(x=vals$ff1, y=vals$sm), 
                type="message")[-c(1:3)], sep="\n"))
        vals$cmp2 <- CATALYST::compCytof(x=vals$ff2, y=vals$sm)
        
    })

    output$plot_plotSpillmat <- renderPlot({
        if (is.null(vals$sm)) return()
        input$button_newSpill
        CATALYST::plotSpillmat(bc_ms=vals$key, SM=isolate(vals$sm))
    })
    
    output$compensation_sidebar_2 <- renderUI({
        if (is.null(vals$cmp1) 
            || is.null(vals$cmp2) 
            || is.null(vals$sm)) 
            return()
        compensation_sidebar_2
    })

# ------------------------------------------------------------------------------
# scatters
# ------------------------------------------------------------------------------
    
    output$text_spill <- renderText({
        input$button_view
        paste0(sprintf("%.3f", 100*vals$sm[
            isolate(input$input_scatterCh1),
            isolate(input$input_scatterCh2)]), "%")
    })
    
    # adjust spill
    observeEvent(input$button_newSpill, {
        if (is.null(input$input_newSpill)) 
            return()
        vals$sm[input$input_scatterCh1, 
                input$input_scatterCh2] <- input$input_newSpill/100
    })

# ------------------------------------------------------------------------------
# download buttons
# ------------------------------------------------------------------------------
    
    output$dwnld_fcs  <- downloadHandler(filename = function()     { "fcs.zip" },
                                         content  = function(file) { tmpdir <- tempdir(); setwd(tmpdir)
                                         inds <- c(0, rownames(bc_key(vals$re2))) %in% sort(unique(bc_ids(vals$re2)))
                                         if (input$box_IDsAsNms) {
                                             CATALYST::outFCS(x=vals$re2, out_path=tmpdir) 
                                             file_nms <- c("Unassigned", rownames(bc_key(vals$re2)))[inds]
                                         } else if (input$box_upldNms) {
                                             if (is.null(input$input_upldNms)) {
                                                 showNotification("Please upload a naming.",
                                                                  type="error", closeButton=FALSE)
                                                 return()
                                             } else if (nrow(vals$nms) < nrow(bc_key(vals$re2))) {
                                                 showNotification(paste("Only", nrow(vals$nms), 
                                                                  "file names provided but", 
                                                                  nrow(bc_key(vals$re2)), "needed."),
                                                                  type="error", closeButton=FALSE)
                                                 return()
                                             } else if (sum(vals$nms[, 1] %in% rownames(bc_key(vals$re2))) != nrow(bc_key(vals$re2))) {
                                                 showNotification("Couldn't find a file name for all samples.\n
                                                                  Please make sure all sample IDs occur\n
                                                                  in the provided naming scheme.",
                                                                  type="error", closeButton=FALSE)
                                                 return()
                                             }
                                             CATALYST::outFCS(x=vals$re2, out_path=tmpdir, out_nms=paste0(vals$nms[, 2], "_", rownames(bc_key(vals$re2))))
                                             file_nms <- c("Unassigned", paste0(vals$nms[, 2], "_", rownames(bc_key(vals$re2))))[inds]
                                         }
                                         zip(zipfile=file, files=paste0(file_nms, ".fcs")) }, 
                                         contentType = "application/zip")
    
    output$dwnld_yep   <- downloadHandler(filename = function()     { "yield_event_plots.zip" },
                                          content  = function(file) { tmpdir <- tempdir(); setwd(tmpdir)
                                          CATALYST::plotYields(x=vals$re2, which=vals$yp_choices, out_path=tmpdir)
                                          CATALYST::plotEvents(x=vals$re2, which=vals$ep_choices, out_path=tmpdir, n_events=as.numeric(input$n_events))
                                          zip(zipfile=file, contentType = "application/zip", files=paste0(c("summary_yield_plot", "yield_plot", "event_plot"), ".pdf")) })
    
    output$dwnld_comped_1 <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs), "-comped.fcs") },
                                             content  = function(file) { flowCore::write.FCS(vals$cmp1, file) })
    output$dwnld_comped_2 <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs2),"-comped.fcs") },
                                             content  = function(file) { flowCore::write.FCS(vals$cmp2, file) })
    output$dwnld_spillMat <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs), "-spillMat.csv") },
                                             content  = function(file) { write.csv(vals$sm, file) })
})
