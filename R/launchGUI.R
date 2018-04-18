#' @rdname launchGUI
#' @title CATALYST Shiny app
#' 
#' @description 
#' Launches the CATALYST Shiny app.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @importFrom shiny runApp
#' @import DT shinyBS shinydashboard shinyjs
#' @export
# ==============================================================================
launchGUI <- function() {
    runApp(
        appDir=system.file("shinyGUI", package="CATALYST"), 
        launch.browser=TRUE)
}