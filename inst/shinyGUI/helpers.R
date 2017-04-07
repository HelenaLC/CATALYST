# ==================================================================================================
# debarcoding tab summary table 
# --------------------------------------------------------------------------------------------------

summary_tbl <- function(x) {
    
    # better way to do this?
    ids <- rownames(x@bc_key)
    yields <- x@yields[cbind(1:nrow(x@bc_key), 
        match(x@sep_cutoffs, seq(0, 1, .01)))]
    tbl <- matrix(c(paste0("<b>", ids, "</b>"), 
        sapply(ids, function(k) sum(x@bc_ids == k)), 
        x@sep_cutoffs, round(yields, 4) * 100), ncol=4)
    colnames(tbl) <- c("IDs", "Counts", "Cutoffs", "Yields")
    cols <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn")[-c(1, 11)])(100)
    if (nrow(x@bc_key) < 21) {
        dom <- "t"
    } else {
        dom <- "tp"
    }
    DT::datatable(tbl, 
                  class="stripe", selection="none", escape=FALSE,
                  extensions="FixedHeader",
                  options=list(scrollY=TRUE, ordering=FALSE, selection=FALSE, 
                               searching=FALSE, info=FALSE, pageLength=20, dom=dom)) %>% 
        DT::formatStyle("Yields", 
                        backgroundColor=DT::styleInterval(cuts=1:99, values=cols))
                            
}

# ==================================================================================================
# scatter plot
# --------------------------------------------------------------------------------------------------

scatter <- function(ff, which, cofactor=50, out_path=NULL, name_ext=NULL) {
    es <- flowCore::exprs(ff)  
    es <- asinh(es / cofactor) 
    if (nrow(ff) > 1e4)         
        es <- es[sample(nrow(ff), 1e4), ]
    
    df <- data.frame(es[, which[1]], es[, which[2]])
    
    lower <- apply(df, 2, function(x)   floor(min(x) * 2) / 2)
    upper <- apply(df, 2, function(x) ceiling(max(x) * 2) / 2)
    
    x_lims <- c(lower[1], upper[1])
    y_lims <- c(lower[2], upper[2])
    
    ggplot(data.frame(es), aes_string(x=which[1], y=which[2])) +
        geom_point(alpha=.1, size=2) + 
        geom_rug(sides="tr", col="darkred", alpha=.2, size=.1) + 
        coord_cartesian(xlim=x_lims, ylim=y_lims, expand=FALSE) +
        scale_x_continuous(labels=function(x) format(x, nsmall=1)) +
        scale_y_continuous(labels=function(x) format(x, nsmall=1)) +
        theme_bw() + theme(
            aspect.ratio=1,
            plot.margin=unit(c(.5,.5,.5,.5), "cm"),
            panel.border=element_rect(colour="black", fill=NA, size=1),
            panel.grid.major=element_line(colour="lightgrey"),
            panel.grid.minor=element_blank(),
            axis.ticks=element_line(size=.5),
            axis.ticks.length=unit(.2, "cm"), 
            axis.title=element_text(size=14), 
            axis.text=element_text(size=12))
}

# ==================================================================================================
# get medians from rectangular brush
# --------------------------------------------------------------------------------------------------

text_info <- function(ff, cofactor, brush, ch1, ch2) {
    if (is.null(brush)) {
        paste("Brush points to\n check medians.")
    } else {
        selected <- brushedPoints(
            df=data.frame(asinh(flowCore::exprs(ff) / cofactor)),
            brush=brush, 
            xvar=ch1,
            yvar=ch2)
        xmed <- sprintf("% .3f", median(selected[, ch1]))
        ymed <- sprintf("% .3f", median(selected[, ch2]))
        paste0(sprintf("%5s", ch1), ": ", xmed, "\n", 
               sprintf("%5s", ch2), ": ", ymed)
    }
}