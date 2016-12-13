# ==============================================================================
# retrieve legend from ggplot
# ------------------------------------------------------------------------------
get_legend <- function(p) {
    tmp <- ggplot_gtable(ggplot_build(p)) 
    lgd <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    return(tmp$grobs[[lgd]]) 
}

# ==============================================================================
# scientific annotation
# ------------------------------------------------------------------------------
scientific_10 <- function(x)
    parse(text=gsub("e", " %*% 10^", format(x, scientific=TRUE)))
