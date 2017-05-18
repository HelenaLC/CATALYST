# load required packages
library(magrittr)

# set maximum web request size to 500 MB
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
    
    source("helpers-normalization_tmp.R")
    library(shinyjs)
    library(shinyBS)
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
        ff_norm = NULL,         # input FCS for normalization
        ffsNorm = list(),       # input flowFrame(s) to be normalized
        ffsNormTo = list(NULL), # flowFrame(s) to use as baseline
        ffsNormed = list(NULL), # normalized flowFrame(s)
        beadGates = list(NULL), # list of bead gates
        beadInds = list(NULL),  # logical vector for each flowFrame
        beads = NULL,    # beads for normalization
        beadCols = NULL, # bead columns
        dnaCols = NULL,  # DNA columns
        timeCol = NULL,  # time column
        smpl = 1,        # index of currently selected sample
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
        ffsComp = list(), # input flowFrame(s) for compensation
        ffsComped = list(NULL),
        log  = NULL, # console output
        trm  = NULL, # trim value
        sm0  = NULL, # original spillover matrix
        sm   = NULL, # spillover matrix
        ch1  = NULL, # scatter channel 1
        ch2  = NULL, # scatter channel 2
        cmp1 = NULL, # compensated flowFrame 1
        cmp2 = NULL) # compensated flowFrame 2


# ==============================================================================
# NORMALIZATION
# ==============================================================================
    
    smpls <- reactive({
        if (input$box_NormToCurrent == 1) {
            input$fcsNorm$name
        } else {
            c(input$fcsNorm$name, input$fcsNormTo$name)
        }
    })
     
    observe({
        if (is.null(vals$beads)) return()
        chs <- flowCore::colnames(ffs()[[vals$smpl]])
        vals$beadCols <- get_bead_cols(chs, isolate(vals$beads))
        vals$dnaCols  <- grep("Ir191|Ir193", chs, ignore.case=TRUE)
        vals$timeCol <- grep("time", chs, ignore.case=TRUE)
    })
    
# ==============================================================================
# DEBARCODING
# ==============================================================================
    
    output$debarcodingSidebar1 <- renderUI({
        if (is.null(input$fcsDeba)) return()
        vals$ff1 <- flowCore::read.FCS(input$fcsDeba$datapath)
        debarcodingSidebar1
    })
    
    # actionButton: Assign preliminary IDs
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
        
        output$debarcodingSidebar2 <- renderUI ({ debarcodingSidebar2 })
    })
    
# ==============================================================================
# COMPENSATION
# ==============================================================================
    
# ------------------------------------------------------------------------------
# checkboxes "Enter trim value" & "Use medians"
# ------------------------------------------------------------------------------
  
    output$text_compCytof <- renderUI({
        x <- input$box_estSM
        if (is.null(x) || x == 0)
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

    # get naming sheet from fileInput
    observe({
        x <- input$input_upldNms
        if (is.null(x) || x == 0) 
            return()
        vals$nms <- read.csv(x$datapath, header=FALSE)
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
    
    output$dwnld_comped_1 <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcsDeba), "-comped.fcs") },
                                             content  = function(file) { flowCore::write.FCS(vals$cmp1, file) })
    output$dwnld_comped_2 <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcs2),"-comped.fcs") },
                                             content  = function(file) { flowCore::write.FCS(vals$cmp2, file) })
    output$dwnld_spillMat <- downloadHandler(filename = function()     { file.path(gsub(".fcs", "", input$fcsDeba), "-spillMat.csv") },
                                             content  = function(file) { write.csv(vals$sm, file) })
})
