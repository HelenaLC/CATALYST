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

# helper to style textInput
textInputCol <- 'shinyjs.textInputCol = function(pars) {
    $("#"+pars.id).css("background-color",pars.col)
    .css("height","24px").css("line-height","12px")
    .css("padding-top","6px").css("padding-bottom","6px")
    .css("margin-bottom","-10px").css("border","0px")
    .css("font-size","12px").css("color","black");}'

# helper to force collapse shinydashboard::box
collapseBox <- "shinyjs.collapse=function(boxId){
    $('#'+boxId).closest('.box').not('.collapsed-box')
    .find('[data-widget=collapse]').click();}"

source("module-yieldPlot.R")
source("module-debaPars.R")

source("ui-normalization.R")   
source("ui-debarcoding.R")
source("ui-compensation.R")