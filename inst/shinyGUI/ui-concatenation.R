# ==============================================================================
# CONCATENATION
# ==============================================================================

concatenationTab <- fluidPage(
    shinyjs::useShinyjs(),
    fluidRow(
        column(
            width=4,
            style="padding:0px", 
            shinydashboard::box(
                title="Upload FCS",
                status="warning",
                solidHeader=TRUE,
                width=12,
                fileInput(
                    inputId="fcsConcat",
                    label=NULL,
                    multiple=TRUE),
                uiOutput("concatOutput")
            )
        ),
        column(
            width=8,
            style="padding:0px", 
            uiOutput("editFcsConcat")
        )
    )
)

concatOutput <- function(outNm) {
    tagList(
        textInput(
            inputId="concatOutputFileNm",
            label="Output file name",
            value=outNm),
        checkboxInput(
            inputId="addFileNo",
            label="Add file number channel"),
        checkboxInput(
            inputId="orderByTime",
            label="Order by time",
            value=TRUE),
        tags$style("#goToNorm {width:49%; float:left}"),
        tags$style("#dwnldConcat {color:white; width:49%; float:right}"),
        div(style="display:inline-block; width:49%; float:left",
            bsButton(
                inputId="goToNorm",
                label="Go to normalization",
                width="100%")),
        downloadButton(
            outputId="dwnldConcat",
            label="Merge files",
            class="btn-success")
    )
}

editFcsConcat <- function(ff) {
    pars <- flowCore::colnames(ff)
    desc <- ff@parameters@data$desc
    n <- ncol(ff)
    parList <- lapply(seq_len(n), function(i) 
        div(style="height:32px; margin-top:0px; margin-bottom:0px", 
            disabled(textInput(
                inputId=paste0("par", i), 
                label=NULL, 
                value=pars[i], 
                width="100%"))))
    desList <- lapply(seq_len(n), function(i) 
        div(style="height:32px; margin-top:0px; margin-bottom:0px", 
            textInput(
                inputId=paste0("des", i),
                label=NULL, 
                value=desc[i], 
                width="100%")))
    shinydashboard::box(
        title="File editing",
        status="warning",
        solidHeader=TRUE,
        width=12,
        wellPanel(
            style="background-color:white; border:0px; 
                overflow-y:scroll; max-height:100vh",
            fluidRow(
                column(
                    width=6,
                    p(strong("Channel name")),
                    do.call(tagList, parList)),
                column(
                    width=6,
                    p(strong("Description")),
                    do.call(tagList, desList)))))
}