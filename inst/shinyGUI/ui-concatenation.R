# ==============================================================================
# CONCATENATION
# ==============================================================================

concatenationTab <- fluidPage(
    fluidRow(
        column(
            width=4,
            shinydashboard::box(
                title="Upload FCS",
                status="warning",
                solidHeader=TRUE,
                width=12,
                fileInput(
                    inputId="fcsConcat",
                    label=NULL
                )
            )
        ),
        column(
            width=8,
            uiOutput("editConcatFcs")
        )
    )
)

editConcatFcs <- shinydashboard::box(
    title="File editing",
    status="warning",
    solidHeader=TRUE,
    width=12
)