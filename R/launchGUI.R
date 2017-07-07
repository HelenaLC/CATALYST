#' @rdname launchGUI
#' @title CATALYST Shiny app
#' 
#' @description 
#' Launches the CATALYST Shiny app.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @import shinydashboard ggplot2
#' @importFrom DT datatable formatStyle
#' @importFrom flowCore read.FCS write.FCS exprs colnames parameters
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly
#' @importFrom shiny runApp
#' @importFrom shinyBS bsButton bsTooltip
#' @importFrom shinyjs useShinyjs extendShinyjs toggle toggleState 
#' @export
# ==============================================================================
launchGUI <- function() {
    runApp(appDir=system.file("shinyGUI", package="CATALYST"), launch.browser=TRUE)
}