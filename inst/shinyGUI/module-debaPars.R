debaParsModule <- function(module) {
    tagList(
        hr(style="border-color:black"),
        # automated cutoff estimation
        # population-specific cutoffs
        # global cutoff
        radioButtons(
            inputId=paste0("deba_cutoffs", module),
            label=NULL,
            choiceValues=c("est_cutoffs", "adj_cutoffs", "global_cutoff"),
            choiceNames=c(
                "Estimate separation cutoffs",
                "Adjust population-specific cutoffs",
                "Enter global separation cutoff")),
        uiOutput(outputId=paste0("deba_cutoffs_UI", module)),
        # mahalanobis distance cutoff
        div(style="display:inline-block; vertical-align:middle; width:75%",
            sliderInput(
                inputId=paste0("mhlCutoff", module),
                label="Mahalanobis distance threshold",
                min=10, 
                max=100, 
                value=30,
                width="100%")),
        div(style="display:inline-block; vertical-align:middle", 
            bsButton(
                inputId=paste0("button_mhlCutoff", module),
                label=NULL,
                icon=icon("share"),
                style="primary",
                size="extra-small")),
        bsTooltip(
            id=paste0("button_mhlCutoff", module),
            title="Apply",
            placement="right"))
}

# ------------------------------------------------------------------------------
# UI for cutoff adjustment
# ------------------------------------------------------------------------------
adjustCutoffUI <- function(dbFrame, choices, module) {
    tagList(
        div(style="display:inline-block; vertical-align:top; width:25%",
            selectInput(
                inputId=paste0("select_adjustCutoff", module),
                label=NULL,
                choices=choices)),
        div(style="display:inline-block; width:25%",
            numericInput(
                inputId=paste0("input_adjustCutoff", module),
                label=NULL,
                value=sep_cutoffs(dbFrame)[1],
                min=0, max=1, step=.01)),
        div(style="display:inline-block; vertical-align:middle",
            tagList(
                bsButton(
                    inputId=paste0("button_adjustCutoff", module),
                    label=NULL,
                    icon=icon("share"),
                    style="primary",
                    size="extra-small"),
                bsTooltip(
                    id=paste0("button_adjustCutoff", module),
                    title="Adjust",
                    placement="right"))))
}

# ------------------------------------------------------------------------------
# UI for global cutoff selection
# ------------------------------------------------------------------------------
globalCutoffUI <- function(module) {
    tagList(
        div(style="display:inline-block; width:25%", 
            numericInput(
                inputId=paste0("input_globalCutoff", module),
                label=NULL, value=NULL, 
                min=0, max=1, step=.01)),
        div(style="display:inline-block; vertical-align:middle",
            bsButton(
                inputId=paste0("button_globalCutoff", module),
                label=NULL, 
                icon=icon("share"), 
                style="primary",
                size="extra-small")),
        bsTooltip(
            id=paste0("button_globalCutoff", module),
            title="Apply", 
            placement="right"))
}