# ==============================================================================
# check validity of input FCS files
# ------------------------------------------------------------------------------

check_FCS_fileInput <- function(input, n=1) {
    invalid <- sum(!sapply(seq_len(n), function(i)
        flowCore::isFCSfile(input[[i, "datapath"]])))
    if (n == 1 && invalid == 1) {
        showNotification(type="error",
            paste0("Input file does not seem to be
                a valid FCS2.0, FCS3.0 or FCS3.1 file."))
        return(FALSE)
    } else if (invalid == 1) {
        showNotification(type="error",
            paste0("1/", n, " files does not seem to be
                a valid FCS2.0, FCS3.0 or FCS3.1 file."))
        return(FALSE)
    } else if (invalid > 1) {
        showNotification(type="error",
            paste0(invalid, "/", n, " files don't seem to be
                valid FCS2.0, FCS3.0 or FCS3.1 files."))
        return(FALSE)
    }
    return(TRUE)
}

# ==============================================================================
# FCS file editing
# ------------------------------------------------------------------------------

alter_pars <- function(ff, pars) {
    es <- exprs(ff)
    colnames(es) <- pars
    ps <- parameters(ff)
    ps@data$name <- pars
    ds <- description(ff)
    ds[paste0("$P", seq_len(ncol(ff)), "N")] <- pars
    flowFrame(es, ps, ds)
}

# ==============================================================================
# get list of masses and metals from measurement parameters
# ------------------------------------------------------------------------------

get_ms_and_mets <- function(chs) {
    ms <- CATALYST:::.get_ms_from_chs(chs)
    mets <- CATALYST:::.get_mets_from_chs(chs)
    setNames(list(ms, mets), c("Mass", "Metal"))
}

# ==============================================================================
# get named list of indices of duplicate masses
# ------------------------------------------------------------------------------

get_duplicate_ms <- function(ms) {
    counts <- table(ms)
    unique <- counts == 1
    if (all(unique)) 
        return()
    d <- as.numeric(names(unique)[!unique])
    d <- d[!is.na(d)]
    n <- length(d)
    setNames(lapply(d, function(i) ms == i), d)
}

# ==============================================================================
# check validity of barcoding scheme CSV
# ------------------------------------------------------------------------------

check_key <- function(key, ff) {
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
    ids <- rownames(bc_key(x))
    yields <- yields(x)[cbind(seq_len(nrow(bc_key(x))), 
        findInterval(sep_cutoffs(x), seq(0, 1, .01)))]
    yields <- sprintf("%2.2f", round(100*yields, 4))
    tbl <- matrix(c(paste0("<b>", ids, "</b>"), sapply(ids, function(k) 
        sum(bc_ids(x) == k)), sep_cutoffs(x), yields), ncol=4,
        dimnames=list(NULL, c("ID", "Count", "Cutoff", "Yield")))
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
