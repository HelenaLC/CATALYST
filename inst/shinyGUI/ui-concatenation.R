# ==============================================================================
# CONCATENATION
# ==============================================================================

concatenationTab <- fluidPage(
    shinyjs::useShinyjs(),
    shinyjs::extendShinyjs(text=restyleTextInputConcat),
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

editFcsConcat <- function(ff) {
    pars <- flowCore::colnames(ff)
    desc <- ff@parameters@data$desc
    n <- ncol(ff)
    parList <- lapply(seq_len(n), function(i)
        textInput(
            inputId=paste0("par", i),
            label=NULL,
            value=pars[i],
            width="100%"))
    desList <- lapply(seq_len(n), function(i)
        textInput(
            inputId=paste0("des", i),
            label=NULL,
            value=desc[i],
            width="100%"))
    shinydashboard::box(
        title="File editing",
        status="warning",
        solidHeader=TRUE,
        width=12,
        wellPanel(
            style="background-color:white; border:0px; 
                overflow-y:scroll; max-height:100vh",
            fluidPage(
            fluidRow(
                column(
                    width=6,
                    p(strong("Channel name")),
                    do.call(tagList, parList)),
                column(
                    width=6,
                    p(strong("Description")),
                    do.call(tagList, desList))
            )
        )
    ))
}