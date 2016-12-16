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

# ==============================================================================
# get spillover columns
# ------------------------------------------------------------------------------
get_spill_cols <- function(ms, mets) {
    
    iso_tbl <- list(
        La=138:139, Pr=141, Nd=c(142:146, 148, 150), 
        Sm=c(144, 147:150, 152, 154), Eu=c(151, 153), 
        Gd=c(152, 154:158, 160), Dy=c(156, 158, 160:164), 
        Tb=159, Er=c(162, 164, 166:168, 170), Ho=165, 
        Yb=c(168, 170:174, 176), Tm=169, Lu=175:176)
    
    spill_cols <- list()
    for (i in seq_along(ms)) {
        if (is.na(ms[i])) next
        p1 <- m1 <- ox <- iso <- NULL
        if ((ms[i] + 1)  %in% ms) p1 <- which(ms == (ms[i] + 1))
        if ((ms[i] - 1)  %in% ms) m1 <- which(ms == (ms[i] - 1)) 
        if ((ms[i] + 16) %in% ms) ox <- which(ms == (ms[i] + 16))
        iso <- iso_tbl[[mets[!is.na(ms)][i]]]
        iso <- which(ms %in% iso[iso != ms[i]])
        spill_cols[[i]] <- unique(c(m1, p1, iso, ox))
    }
    spill_cols
}

# ==============================================================================
# make spillover matrix symmetrical
# ------------------------------------------------------------------------------
make_symetric <- function(SM) {
    sm <- diag(ncol(SM))
    rownames(sm) <- colnames(sm) <- colnames(SM)
    sm[rownames(SM), colnames(SM)] <- SM 
    sm
}
    