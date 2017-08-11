# ==============================================================================
# get spillover columns
# ------------------------------------------------------------------------------
get_spill_cols <- function(ms, mets, l=CATALYST::isotope_list) {
    ms <- as.numeric(ms)
    spill_cols <- vector("list", length(ms))
    for (i in seq_along(ms)) {
        p1 <- m1 <- ox <- iso <- NULL
        if ((ms[i] + 1)  %in% ms) p1 <- which(ms == (ms[i] + 1))
        if ((ms[i] - 1)  %in% ms) m1 <- which(ms == (ms[i] - 1)) 
        if ((ms[i] + 16) %in% ms) ox <- which(ms == (ms[i] + 16))
        iso <- l[[mets[i]]]
        iso <- which(ms %in% iso[iso != ms[i]])
        spill_cols[[i]] <- unique(c(m1, p1, iso, ox))
    }
    spill_cols
}

# ==============================================================================
# compute channel i to j spill
# ------------------------------------------------------------------------------
get_sij <- function(pos_i, neg_i, pos_j, neg_j, method, trim) {
    if (length(neg_i) == 0) neg_i <- 0
    if (length(neg_j) == 0) neg_j <- 0
    if (method == "default") {
        bg_j <- mean(neg_j, trim=.1)
        bg_i <- mean(neg_i, trim=.1)
        receiver <- pos_j - bg_j
        spiller  <- pos_i - bg_i
    } else if (method == "classic") {
        receiver <- mean(pos_j, trim) - mean(neg_j, trim)
        spiller  <- mean(pos_i, trim) - mean(neg_i, trim)
    } else {
        stop("'method = ", method, "' is not a valid option.")
    }
    receiver[receiver < 0] <- 0
    spiller [spiller  < 0] <- 0
    sij <- receiver / spiller
    sij[is.infinite(sij)] <- 0
    median(sij)
}

# ==============================================================================
# make spillover matrix symmetrical
# ------------------------------------------------------------------------------
make_symetric <- function(x) {
    x <- as.matrix(x)
    dims <- dim(x)
    i <- which.max(dim(x))
    nms <- dimnames(x)[[i]]
    M <- matrix(0, dims[i], dims[i])
    rownames(M) <- colnames(M) <- nms
    M[rownames(x), colnames(x)] <- as.matrix(x)
    M
}

# ==============================================================================
# plot for estTrim()
# ------------------------------------------------------------------------------

plot_estTrim <- function(df, trms, xMin, xMax, yMin, yMax, rect, text) {
    opt <- trms[which.min(text$e)]
    ggplot(df, aes_string(x="t", y="m")) +
        geom_vline(aes(xintercept=opt), lty=2, size=.5) +
        geom_jitter(aes_string(fill="Spiller", group="Receiver"),
            col="navy", height=0, width=diff(trms)[1]/5, size=2, alpha=.3) + 
        geom_rect(fill="aliceblue", inherit.aes=FALSE, data=rect, 
            aes_string(xmin="x1", xmax="x2", ymin="y1", ymax="y2")) + 
        geom_text(size=3, col="blue", vjust=.5, data=text, 
            aes_string(label="e", x="x", y="y")) +
        geom_hline(aes(yintercept=0), lty=2, col="red", size=.5) +
        scale_x_continuous(limits=c(xMin, xMax), 
            expand=c(0,0), breaks=trms, labels=format(trms, 2)) +
        scale_y_continuous(limits=c(yMin, yMax+.6), 
            expand=c(0,0), breaks=c(0, yMin:yMax)) +
        labs(x="Trim value used for estimation of spill values", 
            y="Median counts upon compensation") + 
        theme_classic() + 
        theme(legend.position="none", 
            axis.text=element_text(size=8),
            axis.title=element_text(size=10), 
            panel.grid.major.y=element_blank(),
            panel.grid.major.x=element_line(size=.25, color="grey"),
            panel.grid.minor=element_blank())
}

# ==============================================================================
# check validity of input spillover matrix in compCytof()
# ------------------------------------------------------------------------------
check_spillMat <- function(sm) {
    if (any(sm < 0))
        stop("\nThe supplied spillover matrix is invalid ",
            "as it contains negative entries.\n",
            "Valid spill values are non-negative and mustn't exceed 1.")
    if (any(sm > 1))
        stop("\nThe supplied spillover matrix is invalid ",
            "as it contains entries greater than 1.\n",
            "Valid spill values are non-negative and mustn't exceed 1.")
    
    cnames <- colnames(sm)[which(colnames(sm) %in% rownames(sm))]
    sii <- sm[cbind(cnames, cnames)]
    if (any(sii != 1))
        stop("\nThe supplied spillover matrix is invalid ",
            "as its diagonal contains entries != 1.\n")
}


# ==============================================================================
# Helper functions to get mass and metal from a channel name
# ------------------------------------------------------------------------------

get_ms_from_chs <- function(chs) {
    gsub("[[:punct:][:alpha:]]", "", chs)
}

get_mets_from_chs <- function(chs) { 
    gsub("([[:punct:]]*)([[:digit:]]*)(Di)*", "", chs)
}

# ==============================================================================
# This function compares a list of spillover channels to
# a provided spillover matrix and provides warnings for cases
# where the spillovermatrix has no information about a potential
# expected spillover among the new channels.
# ------------------------------------------------------------------------------

warn_new_intearctions <- function(chs_new, sm) {
    chs_emitting  <- rownames(sm)
    chs_receiving <- colnames(sm)
    chs <- list(chs_new, chs_emitting, chs_receiving)
    chs <- setNames(chs, c("new", "emitting", "receiving"))
    
    # get the metals and masses from the names
    mets <- lapply(chs, get_mets_from_chs)
    ms <- lapply(chs, get_ms_from_chs)
    
    # get the potential mass channels a channel could cause spillover in
    spill_cols <- get_spill_cols(ms$new, mets$new)
    
    first <- TRUE
    for (i in order(ms$new)) {
        # check if the provided spillovermatrix had the 
        # current channel measured as a spill emitter
        is_new_emitting <- !chs$new[i] %in% chs$emitting
        cur_spillms <- ms$new[spill_cols[[i]]]
        
        if (is_new_emitting) {
            # if channel was never measured, no information is present
            # for spill in any of the potential spill receiver channels
            mass_new_rec <- cur_spillms
        } else {
            # if it has been measured, no information is only present
            # for the channels which were not considered spill receivers 
            # in the spillover matrix provided
            mass_new_rec <- cur_spillms[!cur_spillms %in% ms$receiving]
        }
        
        if (length(mass_new_rec) > 0) {
            if (first) {
                message("WARNING: ",
                    "Compensation is likely to be inaccurate.\n",
                    "         ",
                    "Spill values for the following interactions\n",
                    "         ",
                    "have not been estimated:")
                first <- FALSE
            }
            message(chs$new[i], " -> ", paste(
                chs$new[ms$new %in% mass_new_rec], collapse=", "))
        }
    }
}