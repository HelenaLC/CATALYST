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

# guides 
    source("ui/guides.R")
    output$debarcoding_guide  <- renderUI(debarcoding_guide)
    output$compensation_guide <- renderUI(compensation_guide)
    
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
        beads = NULL,    # beads for normalization
        beadInds = list(NULL),
        smpl = 1,        # index of currently selected sample
# ----- DEBARCODING -----
        ffsDeba = NULL, # input flowFrame(s) for debarcoding
        DebaKey = NULL, # debarcoding scheme
        reAssignPrelim = NULL, # result assignPrelim()
        reEstCutoffs = NULL,   # result estCutoffs()
        reApplyCutoffs = NULL, # result applyCutoffs()
        mhlCutoff = 30,     # mhl_cutoff for applyCutoffs()
        ids = NULL,         # current IDs with event assignments
        yp_choices  = NULL, # all barcode IDs including 0
        ep_choices  = NULL, # only IDs with events assigned
        adj_choices = NULL, # all barcode IDs excluding 0
        mhl_choices = NULL, # IDs with events assigned excluding 0
        nms = NULL,         # names to use for output FCS files
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
    
    output$dwnld_fcs <- downloadHandler(
        filename=function() { "fcs.zip" },
        content =function(file) { 
            tmpdir <- tempdir()
            setwd(tmpdir)
            inds <- c(0, rownames(bc_key(vals$reApplyCutoffs))) %in% 
                sort(unique(bc_ids(vals$reApplyCutoffs)))
            if (input$box_IDsAsNms) {
                CATALYST::outFCS(x=vals$reApplyCutoffs, out_path=tmpdir) 
                file_nms <- c("Unassigned", rownames(bc_key(vals$reApplyCutoffs)))[inds]
            } else if (input$box_upldNms) {
                if (is.null(input$input_upldNms)) {
                    showNotification("Please upload a naming.",
                        type="error", closeButton=FALSE)
                    return()
                } else if (nrow(vals$nms) < nrow(bc_key(vals$reApplyCutoffs))) {
                    showNotification(paste("Only", nrow(vals$nms), 
                        "file names provided but", 
                        nrow(bc_key(vals$reApplyCutoffs)), "needed."),
                        type="error", closeButton=FALSE)
                    return()
                } else if (sum(vals$nms[, 1] %in% rownames(bc_key(vals$reApplyCutoffs))) != nrow(bc_key(vals$reApplyCutoffs))) {
                    showNotification(
                        "Couldn't find a file name for all samples.\n
                        Please make sure all sample IDs occur\n
                        in the provided naming scheme.",
                        type="error", closeButton=FALSE)
                    return()
                }
                CATALYST::outFCS(
                    x=vals$reApplyCutoffs, 
                    out_path=tmpdir, 
                    out_nms=paste0(vals$nms[, 2], "_", rownames(bc_key(vals$reApplyCutoffs))))
                file_nms <- c("Unassigned", paste0(vals$nms[, 2], "_", rownames(bc_key(vals$reApplyCutoffs))))[inds]
            }
            zip(zipfile=file, files=paste0(file_nms, ".fcs")) }, 
        contentType = "application/zip")                        
    
    output$dwnld_yep <- downloadHandler(
        filename=function() { "yield_event_plots.zip" },
        content =function(file) { 
            tmpdir <- tempdir()
            setwd(tmpdir)
            CATALYST::plotYields(
                x=vals$reApplyCutoffs, 
                which=vals$yp_choices, 
                out_path=tmpdir)
            CATALYST::plotEvents(
                x=vals$reApplyCutoffs, 
                which=vals$ep_choices, 
                out_path=tmpdir, 
                n_events=as.numeric(input$n_events))
            zip(zipfile=file, 
                contentType="application/zip", 
                files=paste0(c("summary_yield_plot", "yield_plot", "event_plot"), ".pdf")) 
        }
    )
    
    output$dwnld_comped_1 <- downloadHandler(
        filename=function() { 
            file.path(gsub(".fcs", "", input$fcsDeba, ignore.case=TRUE), "_comped.fcs")
        },
        content =function(file) { 
            flowCore::write.FCS(vals$cmp1, file) 
        }
    )
    output$dwnld_comped_2 <- downloadHandler(
        filename=function() { 
            file.path(gsub(".fcs", "", input$fcs2, ignore.case=TRUE), "_comped.fcs") 
        },
        content =function(file) { 
            flowCore::write.FCS(vals$cmp2, file) 
        }
    )
    output$dwnld_spillMat <- downloadHandler(
        filename=function() { 
            file.path(gsub(".fcs", "", input$fcsDeba), "_spillMat.csv") 
        },
        content=function(file) { 
            write.csv(vals$sm, file) 
        }
    )
})
