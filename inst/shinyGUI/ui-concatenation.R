# ==============================================================================
# CONCATENATION
# ==============================================================================

concatenationTab <- fluidPage(
    fluidRow(
        column(
            width=3,
            shinydashboard::box(
                title="Upload FCS",
                status="warning",
                solidHeader=TRUE,
                width=NULL,
                fileInput(
                    inputId="fcsConcat",
                    label=NULL
                )
            )
        )
    )
)