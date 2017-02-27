library(shinydashboard)

header <- dashboardHeader(title=tagList("CATALYST", icon("rocket")))

header$children[[3]]$children <- tagList(
    tags$a(href="http://www.imls.uzh.ch/en/research/robinson.html",
           img(src="http://www.isigrowth.eu/wp-content/uploads/2015/09/LOGO_Zurich2.png", 
               align="right", height="45px", style="padding-top:5px; padding-left:5px; padding-right:30px;")))

source("ui/aes.R")
source("ui/debarcoding_tab.R")
source("ui/compensation_tab.R")

sidebar <- dashboardSidebar(sidebarMenu(
    menuItem("Debarcoding",  tabName="debarcoding",  icon=icon("barcode")),
    menuItem("Compensation", tabName="compensation", icon=icon("cogs"))))

body <- dashboardBody(tabItems(
    tabItem(tabName="debarcoding",  debarcoding_tab),
    tabItem(tabName="compensation", compensation_tab)))
    
dashboardPage(header, sidebar, body, skin="black")