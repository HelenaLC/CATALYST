# load required packages^
library(magrittr)

# set maximum web request size to 500 MB
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {

    source("helpers-normalization_tmp.R")
    library(shinyjs)
    library(shinyBS)
    library(plotly)
    library(flowCore)
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
        customBeads = NULL,
        ff_norm = NULL,         # input FCS for normalization
        beads = NULL,    # beads for normalization
        beadInds = list(NULL),
        smpl = 1,        # index of currently selected sample
# debarcoding
        dbFrame1 = NULL, # preliminary dbFrame
        dbFrame2 = NULL, # dbFrame with deconvolution parameters applied
# compensation
        sm  = NULL,  # original spillover matrix
        cmp1 = NULL, # compensated flowFrame 1
        cmp2 = NULL) # compensated flowFrame 2

# ==============================================================================
# COMPENSATION
# ==============================================================================
    
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
    
    output$dwnld_debaFcs <- downloadHandler(
        filename=function() { "fcs.zip" },
        content =function(file) { 
            ids <- rownames(bc_key(dbFrame()))
            tmpdir <- tempdir()
            setwd(tmpdir)
            inds <- c(0, ids) %in% sort(unique(bc_ids(dbFrame())))
            if (input$box_IDsAsNms) {
                CATALYST::outFCS(x=dbFrame(), out_path=tmpdir) 
                fileNms <- c("Unassigned", ids)[inds]
            } else if (input$box_upldNms) {
                if (is.null(input$input_upldNms)) {
                    showNotification("Please upload a naming.",
                        type="error", closeButton=FALSE)
                    return()
                } else if (nrow(vals$nms) < nrow(bc_key(dbFrame()))) {
                    showNotification(paste("Only", nrow(vals$nms), 
                        "file names provided but", 
                        nrow(bc_key(dbFrame())), "needed."),
                        type="error", closeButton=FALSE)
                    return()
                } else if (sum(vals$nms[, 1] %in% ids) != nrow(bc_key(dbFrame()))) {
                    showNotification(
                        "Couldn't find a file name for all samples.\n
                        Please make sure all sample IDs occur\n
                        in the provided naming scheme.",
                        type="error", closeButton=FALSE)
                    return()
                }
                CATALYST::outFCS(
                    x=dbFrame(), 
                    out_path=tmpdir, 
                    out_nms=paste0(vals$nms[, 2], "_", ids))
                fileNms <- c("Unassigned", paste0(vals$nms[, 2], "_", ids))[inds]
            }
            zip(zipfile=file, files=paste0(fileNms, ".fcs")) }, 
        contentType = "application/zip")                        
    
    output$dwnld_debaPlots <- downloadHandler(
        filename=function() { "yield_event_plots.zip" },
        content =function(file) { 
            tmpdir <- tempdir()
            setwd(tmpdir)
            CATALYST::plotYields(
                x=dbFrame(), 
                which=vals$yp_choices, 
                out_path=tmpdir)
            CATALYST::plotEvents(
                x=dbFrame(), 
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
