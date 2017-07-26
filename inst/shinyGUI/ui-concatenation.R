# ==============================================================================
# CONCATENATION
# ==============================================================================

concatenationTab <- fluidPage(
    fluidRow(
        column(
            width=3,
            shinydashboard::box(
                title="Upload FCS",
                fileInput(
                    inputId="fcsConcat",
                    label=NULL
                )
            )
        )
    )
)