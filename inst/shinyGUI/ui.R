header <- dashboardHeader(title=tagList("CATALYST", icon("rocket")))

header$children[[3]]$children <- tagList(
    tags$a(href="http://www.imls.uzh.ch/en/research/robinson.html",
           img(src="http://www.isigrowth.eu/wp-content/uploads/2015/09/LOGO_Zurich2.png", 
               style="padding-top:5px; padding-left:5px; padding-right:30px;",
               align="right", height="45px")))

sidebar <- dashboardSidebar(sidebarMenu(
    menuItem("Normalization", tabName="normalization", icon=icon("balance-scale")),
    menuItem("Debarcoding",   tabName="debarcoding",   icon=icon("barcode")),
    menuItem("Compensation",  tabName="compensation",  icon=icon("cogs"))))

body <- dashboardBody(tabItems(
    tabItem(tabName="normalization", normalization_tab),
    tabItem(tabName="debarcoding",   debarcoding_tab),
    tabItem(tabName="compensation",  compensation_tab)))
    
dashboardPage(header, sidebar, body, skin="black")