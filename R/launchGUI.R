#' @export
launchGUI <- function() {
    shiny::runApp(appDir=file.path(system.file(package="CATALYST"), "shinyGUI"), display.mode="normal")
}