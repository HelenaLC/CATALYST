editMs <- function(id) {
    # id := 'SS' or 'MP'
    modalDialog(
        fluidPage(
            column(
                offset=2,
                width=4,
                align="center",
                uiOutput(paste0("duplicateMasses", id))),
            column(
                width=4,
                align="center",
                uiOutput(paste0("duplicateMetals", id)))),
        title=HTML(paste(strong("Duplicate masses detected"), 
            "Please select which channels to keep.", sep="<br/>")),
        footer=shinyBS::bsButton(
            inputId=paste0("msChecked", id), 
            label="Done"),
        size="m")
}

editMets <- function(id, header1, header2) {
    # id := 'FLvsSS' or 'SSvsMP'
    # header1 := 'Fulidigm' or 'Single-stains'
    # header2 := 'Single-stains' or 'Multiplexed'
    modalDialog(
        wellPanel(
            style="background-color:white; border:0px; 
            overflow-y:scroll; max-height:100vh",
            fluidPage(
                fluidRow(
                    column(
                        width=4,
                        p(strong("Mass"))),
                    column(
                        width=4,
                        p(strong(header1))),
                    column(
                        width=4,
                        p(strong(header2)))),
                fluidRow(
                    column(
                        width=4,
                        align="center",
                        uiOutput(paste0("masses", id))),
                    column(
                        width=3,
                        align="center",
                        uiOutput(paste0("metals1", id))),
                    column(
                        width=1,
                        uiOutput(paste0("boxes1", id))),
                    column(
                        width=3,
                        align="center",
                        uiOutput(paste0("metals2", id))),
                    column(
                        width=1,
                        uiOutput(paste0("boxes2", id)))))),
        title=HTML(paste(sep="<br/>", strong("Metal mismatches detected"), 
            "The metals selected here will effect spillover estimation.
            Please check carefully and comfirm which metals were used.")),
        footer=shinyBS::bsButton(
            inputId=paste0("metsChecked", id), 
            label="Done"),
        size="m")
}

msList <- function(n, id) {
    # n := number of channels
    # id := 'FLvsSS' or 'SSvsMP'
    l <- lapply(seq_len(n), function(i) {
        div(style="height:40px; width:100%", 
            verbatimTextOutput(
                output=paste0("mass", id, i)))
    })
    do.call(tagList, l)
}

selectMs <- function(mets, inds, id) {
    l <- lapply(seq_along(inds), function(i) {
        selectInput(
            inputId=paste0("metalSelection", id, i),
            label=NULL,
            choices=mets[inds[[i]]],
            width="100%")
    })
    do.call(tagList, l)
}

metsList <- function(n, id1, id2) {
    # n := number of channels
    # id1 := 'FLvsSS' or 'SSvsMP'
    # id2 := 'ref' or 'dat'
    l <- lapply(seq_len(n), function(i) {
        div(style="height:40px; width:100%", 
            verbatimTextOutput(
                output=paste0("metal", id1, id2, i)))
    })
    do.call(tagList, l)
}

boxList <- function(n, id1, id2) {
    # n := number of channels
    # id1 := 'FLvsSS' or 'SSvsMP'
    # id2 := 'ref' or 'dat'
    val <- id2 == "dat"
    l <- lapply(seq_len(n), function(i) {
        div(style="height:30px; margin:0; padding:0",
            checkboxInput(
                inputId=paste0("box", id1, id2, i),
                label=NULL,
                value=val))
    })
    do.call(tagList, l)
}
