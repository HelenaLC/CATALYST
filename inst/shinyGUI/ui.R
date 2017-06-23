header <- dashboardHeader(title=tagList("CATALYST", icon("rocket")))

header$children[[3]]$children <- tagList(
    tags$a(href="http://www.imls.uzh.ch/en/research/robinson.html",
           img(src="http://www.isigrowth.eu/wp-content/uploads/2015/09/LOGO_Zurich2.png", 
               style="padding-top:5px; padding-left:5px; padding-right:30px;",
               align="right", 
               height="45px")
    )
)

sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem(
            text="Normalization", 
            tabName="normalization", 
            icon=icon("balance-scale")),
        menuSubItem(
            text="Bead gating", 
            tabName=NULL, 
            icon=icon("angle-double-right")),
        menuSubItem(
            text="Bead removal",
            tabName=NULL,
            icon=icon("angle-double-right")),
        menuItem(
            text="Debarcoding",   
            tabName="debarcoding",   
            icon=icon("barcode")),
        menuItem(
            text="Compensation",  
            tabName="compensation",  
            icon=icon("cogs"))
    )
)

body <- dashboardBody(
    tabItems(
        tabItem(tabName="normalization", normalizationTab),
        tabItem(tabName="debarcoding",   debarcodingTab),
        tabItem(tabName="compensation",  compensationTab)
    )
)
    
dashboardPage(header, sidebar, body, skin="black")