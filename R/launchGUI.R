#' @rdname launchGUI
#' @title Launch shiny GUI
#' @author Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' 
launchGUI <- function() {
    runApp(appDir="https://github.com/HelenaLC/CATALYST/shinyGUI", launch.browser=TRUE)
}
