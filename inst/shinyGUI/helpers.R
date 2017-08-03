# ==============================================================================
# check validity of input FCS files
# ------------------------------------------------------------------------------

check_FCS_fileInput <- function(input, n) {
    invalid <- sum(!sapply(seq_len(n), function(i)
        flowCore::isFCSfile(input[[i, "datapath"]])))
    if (n == 1 && invalid == 1) {
        showNotification(
            paste0("Input file does not seem to be
                a valid FCS2.0, FCS3.0 or FCS3.1 file."),
            duration=NULL, type="error")
        return(FALSE)
    } else if (invalid == 1) {
        showNotification(
            paste0("1/", n, " files does not seem to be
                a valid FCS2.0, FCS3.0 or FCS3.1 file."),
            duration=NULL, type="error")
        return(FALSE)
    } else if (invalid > 1) {
        showNotification(
            paste0(invalid, "/", n, " files don't seem to be
                valid FCS2.0, FCS3.0 or FCS3.1 files."),
            duration=NULL, type="error")
        return(FALSE)
    }
    return(TRUE)
}

# ==============================================================================
# FCS file editing
# ------------------------------------------------------------------------------

new_ff <- function(data, pars, desc) {
    colnames(data) <- pars
    params <- data.frame(list(name=pars, desc=desc))
    params$minRange <- apply(data, 2, function(x) min(min(x), 0))
    params$maxRange <- apply(data, 2, max)
    params$range <- params$maxRange - params$minRange
    params <- Biobase::AnnotatedDataFrame(params)
    row.names(params) <- paste0('$P', seq_len(nrow(params)))

    desc <- list()
    for (i in seq_len(nrow(params))){
        nm <- row.names(params)[i]
        desc[paste0(nm, 'N')] <- as.character(params$name[i])
        desc[paste0(nm, 'S')] <- as.character(params$desc[i])
        desc[paste0(nm, 'R')] <- params$range[i]
    }
    flowCore::flowFrame(data, params, desc)
}

# ==============================================================================
# check validity of barcoding scheme CSV
# ------------------------------------------------------------------------------

checkKey <- function(key, ff) {
    if (any(rownames(key) == "") || 
            length(unique(rownames(key))) != nrow(key)) {
        showNotification(
            "BARCODING SCHEME INVALID:
            All samples need to be given a 
            non-empty and unique name.",
            type="error", duration=NULL)
        return(FALSE)
    }
    test <- tryCatch(CATALYST::assignPrelim(
        x=ff[1:10, ], y=key, verbose=FALSE), error=function(e) e)
    if (inherits(test, "error")) {
        showNotification(
            "BARCODING SCHEME INVALID:
            Column names need to be numeric masses
            and match with measurement parameters.",
            type="error", duration=NULL)
        return(FALSE)
    } 
    if (any(!unlist(c(key)) %in% c(0, 1))) {
        showNotification(
            "BARCODING SCHEME INVALID:
            Only binary values are allowed.",
            type="error", duration=NULL)
        return(FALSE)
    }
    return(TRUE)
}

# ==============================================================================
# debarcoding summary table: 
# IDs | Counts | Cutoffs | Yields
# ------------------------------------------------------------------------------

summary_tbl <- function(x) {
    
    # better way to do this?
    ids <- rownames(bc_key(x))
    yields <- yields(x)[cbind(seq_len(nrow(bc_key(x))), 
        findInterval(sep_cutoffs(x), seq(0, 1, .01)))]
    tbl <- matrix(c(paste0("<b>", ids, "</b>"), 
        sapply(ids, function(k) sum(bc_ids(x) == k)), 
        sep_cutoffs(x), round(yields, 4) * 100), ncol=4)
    colnames(tbl) <- c("ID", "Count", "Cutoff", "Yield")
    cols <- colorRampPalette(
        RColorBrewer::brewer.pal(11, "RdYlGn")[-c(1, 11)])(100)
    if (nrow(bc_key(x)) < 21) {
        dom <- "t"
    } else {
        dom <- "tp"
    }
    DT::datatable(tbl, class="stripe", selection="none", 
        escape=FALSE, extensions="FixedHeader",
        options=list(scrollY=TRUE, ordering=FALSE, selection=FALSE, 
            searching=FALSE, info=FALSE, pageLength=20, dom=dom)) %>% 
        DT::formatStyle("Yield", 
            backgroundColor=DT::styleInterval(cuts=1:99, values=cols))
}

# ==============================================================================
# get medians from rectangular brush
# ------------------------------------------------------------------------------

text_info <- function(ff, cf, brush, ch1, ch2) {
    if (is.null(brush)) {
        paste("Brush points to\ncheck medians.")
    } else {
        selected <- brushedPoints(
            df=data.frame(asinh(flowCore::exprs(ff)/cf)),
            brush=brush, 
            xvar=ch1,
            yvar=ch2)
        xmed <- sprintf("% .3f", median(sinh(selected[, ch1])*cf))
        ymed <- sprintf("% .3f", median(sinh(selected[, ch2])*cf))
        paste0(sprintf("%5s", ch1), ": ", xmed, "\n", 
               sprintf("%5s", ch2), ": ", ymed)
    }
}