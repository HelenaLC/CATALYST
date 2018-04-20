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
#' \dontrun{launchGUI}
#'
#' @import shinydashboard 
#' @importFrom shiny runApp
#' @importFrom DT datatable dataTableOutput formatStyle renderDataTable
#' @importFrom shinyBS bsButton bsPopover
#' @importFrom shinyjs hide disable extendShinyjs useShinyjs
#' @export

launchGUI <- function() {
    runApp(
        appDir=system.file("shinyGUI", package="CATALYST"), 
        launch.browser=TRUE)
}