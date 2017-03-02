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
    
   

