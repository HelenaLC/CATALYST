library(CATALYST)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(ggplot2)
library(grid)
library(gridExtra)
library(plotly)
library(flowCore)
library(magrittr)
library(DT)

# helper to force collapse shinydashboard::box
collapseBox <- "shinyjs.collapse=function(boxId){
$('#'+boxId).closest('.box').not('.collapsed-box')
.find('[data-widget=collapse]').click();}"

source("module-yieldPlot.R")
source("module-debaPars.R")

source("ui-normalization.R")   
source("ui-debarcoding.R")
source("ui-compensation.R")