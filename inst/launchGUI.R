#' @import shiny shinydashboard magrittr ggplot2
#' @importFrom DT datatable formatStyle
#' @export
launchGUI <- function() {
    runApp(appDir=file.path(system.file(package="CATALYST"), "shinyGUI"), launch.browser=FALSE)
}