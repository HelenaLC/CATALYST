#' @rdname launchGUI
#' @title Launch GUI
#' 
#' @description Launches the CATALYST Shiny app.
#' 
#' @details Detailed user guides are available inside the app. To use the app 
#' online, please visit \emph{http://imlspenticton.uzh.ch:3838/CATALYSTLite}.
#' 
#' @return Opens a browser window with an interactive Shiny application.
#' 
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
#' @examples
#' launchGUI()
#'
#' @importFrom shiny runApp
#' @import DT shinyBS shinydashboard shinyjs
#' @export

launchGUI <- function() {
    runApp(
        appDir=system.file("shinyGUI", package="CATALYST"), 
        launch.browser=TRUE)
}