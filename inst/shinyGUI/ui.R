header <- dashboardHeader(
    title=tagList("CATALYST", icon("rocket"))

)

header$children[[3]]$children <- tagList(
    tags$a(
        href="http://www.imls.uzh.ch/en/research/robinson.html",
        img(src="http://www.isigrowth.eu/wp-content/uploads/2015/09/LOGO_Zurich2.png", 
            style="padding-top:5px; padding-bottom:5px; padding-left:15px; padding-right:30px;",
            align="right", 
            height="50px")
    ),
    tags$a(
        href="http://bioconductor.org/packages/CATALYST/",
        img(src="https://www.bioconductor.org/images/logo/jpg/bioconductor_logo_rgb.jpg", 
            style="padding-top:5px; padding-bottom:5px; padding-right:15px;",
            background="blue",
            align="right", 
            height="50px")
    ),
    tags$a(
        href="https://github.com/Bioconductor-mirror/CATALYST/tree/release-3.5",
        img(src="https://assets-cdn.github.com/images/modules/logos_page/GitHub-Mark.png",
            style="padding-right:20px;",
            align="right", 
            height="50px")
    )
)

sidebar <- dashboardSidebar(
    sidebarMenu(
        id="tabs",
        menuItem(
            text="Guides",
            tabName="guides",
            icon=icon("question-circle")),
        menuItem(
            text="Concatenation",
            tabName="concatenation",
            icon=icon("clone")),
        menuItem(
            text="Normalization", 
            tabName="normalization", 
            icon=icon("balance-scale")),
        menuItem(
            text="Compensation",  
            tabName="compensation",  
            icon=icon("cogs")),
        menuItem(
            text="Debarcoding",   
            tabName="debarcoding",   
            icon=icon("barcode"))
    )
)

body <- dashboardBody(
    tabItems(
        tabItem(tabName="guides",        guidesTab),
        tabItem(tabName="concatenation", concatenationTab),
        tabItem(tabName="normalization", normalizationTab),
        tabItem(tabName="compensation",  compensationTab),
        tabItem(tabName="debarcoding",   debarcodingTab)
    )
)
    
dashboardPage(header, sidebar, body, skin="black")
