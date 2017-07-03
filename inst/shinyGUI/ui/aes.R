inline       <- "display:inline-block;"
inlineTop    <- "display:inline-block; vertical-align:top"
inlineTop34  <- "display:inline-block; vertical-align:top; height:34px"
inlineCenter <- "display:inline-block; vertical-align:center"

collapseBox <- "shinyjs.collapse=function(boxId){
    $('#'+boxId).closest('.box').find('[data-widget=collapse]').click();}"