# ==============================================================================
# get spillover columns
# ------------------------------------------------------------------------------
get_spill_cols <- function(ms, mets, l=CATALYST::isotope_list) {
    ms <- as.numeric(ms)
    spill_cols <- list()
    for (i in seq_along(ms)) {
        if (is.na(ms[i])) next
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
        theme_classic() + theme(legend.position="none", 
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
get_mass_from_channel <- function(channel){
    return( as.numeric(gsub("[[:punct:][:alpha:]]", "", channel)))
}

get_metal_from_channel <- function(channel){
    return(gsub("[[:digit:]].*", "", channel))
}

# ==============================================================================
# This function compares a list of spillover channels to
# a provided spillover matrix and provides warnings for cases
# where the spillovermatrix has no information about a potential
# expected spillover among the new channels.
# ------------------------------------------------------------------------------

warn_new_intearctions <- function(chan_new, sm){
    chan_emmiting <- rownames(sm)
    chan_receiving <- colnames(sm)
    channels <- list(chan_new, chan_emmiting, chan_receiving)
    
    # get the metals from the names
    metals <- lapply(channels, get_metal_from_channel)
    metals_new <- metals[[1]]
    metals_emmiting <- metals[[2]]
    metals_receiving <- metals[[3]]
    
    # get the mass from the names
    mass = lapply(channels, get_mass_from_channel)
    mass_new <- mass[[1]]
    mass_emmiting <- mass[[2]]
    mass_receiving <- mass[[3]]
    
    # get the potential mass channels a channel could cause spillover in
    spill_cols <- get_spill_cols(mass_new, metals_new)
    
    first <- TRUE
    
    for (i in order(mass_new)) {
        # check if the provided spillovermatrix had the 
        # current channel measured as a spill emitter
        is_new_emmitting <-(!chan_new[i] %in% chan_emmiting)
        cur_spillmasses <- mass_new[spill_cols[[i]]]
        
        if (is_new_emmitting){
            # if channel was never measured, no information is present
            # for spill in any of the potential spill receiver channels
            mass_new_rec = cur_spillmasses
        } else {
            # if it has been measured, no information is only present
            # for the channels which were not considered spill receivers 
            # in the spillover matrix provided
            mass_new_rec <- cur_spillmasses[!cur_spillmasses%in% mass_receiving]
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
            message(chan_new[i], " -> ", paste(
                chan_new[mass_new %in% mass_new_rec], collapse=", "))
        }
        
    }
}

# ==============================================================================
# check which channels of input flowFrame are not 
# contained in spillover matrix and give warning
# ------------------------------------------------------------------------------

prep_spillMat <- function(new_chs, sm_chs) {
    new_mets <- gsub("[[:digit:]]+Di", "", new_chs)
    old_ms <- as.numeric(gsub("[[:punct:][:alpha:]]", "", sm_chs))
    new_ms <- as.numeric(gsub("[[:punct:][:alpha:]]", "", new_chs))
    ms <- c(old_ms, new_ms)
    o <- order(ms)
    ms <- ms[o]
    nms <- c(sm_chs, new_chs)[o]
    # get potential spillover interactions 
    all_mets <- gsub("[[:digit:]]+Di", "", nms)
    spill_cols <- get_spill_cols(ms, all_mets)
    
    first <- TRUE
    for (i in seq_along(new_ms)) {
        idx <- which(ms == new_ms[i] & all_mets == new_mets[i])
        if (length(idx) > 0) {
            if (first) {
                message("WARNING: ",
                        "Compensation is likely to be inaccurate.\n",
                        "         ",
                        "Spill values for the following interactions\n",
                        "         ",
                        "have not been estimated:")
                first <- FALSE
            }
            message(nms[idx], " -> ", paste(
                nms[spill_cols[[idx]]], collapse=", "))
        }
    }
    return(list(old=old_ms, new=new_ms))
}


# ==============================================================================
# Adapt spillover matrix for compensation
# ------------------------------------------------------------------------------

#' @rdname adaptCompensationSpillmat
#' @title Adapts a spillovermatrix such that it can be used directly for compensation
#' 
#' @description 
#' This helper function adapts the columns of a provided spillover matrix such that it is compatible with
#' data having the column names provided.
#'
#' @param input_sm     
#' A previously calculated spillover matrix.
#' @param output_names
#' The names that the prepared output spillover matrix should have as column names.
#' Numeric names as well as names of the form MetalMassDi (e.g. Ir191Di or Ir191) will be
#' interpreted as masses with associated metals.
#' @details
#' The rules how the spillover matrix is adapted can be found in the compCytof function documentation.
#' 
#' @return 
#' An adapted spillovermatrix with column and row names according to the output names
#' 
#' @examples
#' data(ss_exp)
#' 
#' # specify mass channels stained for
#' bc_ms <- c(139, 141:156, 158:176)
#'  re <- assignPrelim(x = ss_exp, y = bc_ms)
#' re <- estCutoffs(x = re)
#' re <- applyCutoffs(x = re)
#' spillMat <- computeSpillmat(x = re)
#' ff <- 
#' adaptCompensationSpillmat(spillMat, colnames(ss_exp))
#'
#' @author 
#' Helena Lucia Crowell \email{crowellh@student.ethz.ch}
#' and Vito Zanotelli \email{vito.zanotelli@uzh.ch}
#' @importFrom flowCore flowFrame colnames exprs compensate
#' @export
# ------------------------------------------------------------------------------

adaptCompensationSpillmat <- function(input_sm, out_channels){
    
    # check the spillovermatrix for obvious errors
    check_spillMat(input_sm)
    # make it symmetric
    #input_sm <- make_symetric(input_sm)
    # get the output names, metals and masses
    n = length(out_channels)
    out_masses <- get_mass_from_channel(out_channels)
    out_metalchannels <- out_channels[!is.na(out_masses)]
    input_sm_channels_col <- colnames(input_sm)
    input_sm_channels_row <- rownames(input_sm)
    input_sm_channels <- unique(c(input_sm_channels_row,
                                  input_sm_channels_col))
    
    # copy the existing spillover information into a new spillover matrix
    sm <- matrix(diag(n), n, n, dimnames=list(out_channels, out_channels))
    sm_preexisting_col <- input_sm_channels_col[input_sm_channels_col %in% out_channels]
    sm_preexisting_row <- input_sm_channels_row[input_sm_channels_row %in% out_channels]
    sm[sm_preexisting_row, sm_preexisting_col] <- input_sm[sm_preexisting_row, sm_preexisting_col]
    
    warn_new_intearctions(out_metalchannels, input_sm)
    
    # check for new channels
    new_metalchannels <- out_metalchannels[!out_metalchannels %in% input_sm_channels]
    new_masses <- get_mass_from_channel(new_metalchannels)
    
    old_receiving_masses <- get_mass_from_channel(input_sm_channels_col)
    
    test <- (length(new_metalchannels) != 0) && (any(inds <- old_receiving_masses %in% new_masses))
    if (test) {
        # check if any new masses were already present in the old masses
        # and add them to receive spillover according to the old masses
        
        # get the channels that correspond to the old_masses 
        # that have an aditional metal with the same weight
        y_col <- input_sm_channels_col[inds]
        names(y_col) <- as.character(old_receiving_masses[inds])
        # get all columns that are part of the affected masses
        fil <- out_masses %in% old_receiving_masses[inds]
        sm_col <- out_channels[fil]
        sm_col_ms <- as.character(out_masses[fil])
        # add the spillover
        old_rowchannels = out_channels[out_channels %in% input_sm_channels_row]
        sm[old_rowchannels, sm_col] <- input_sm[old_rowchannels, y_col[sm_col_ms]]
        for (m in unique(sm_col_ms)){
            mfil <- out_masses == m
            # set the spillover between channels of the same mass to 0
            # otherwise the linear system can get singular
            # diagonal elements will be set to 1 again later on
            sm[mfil, mfil] <- 0
        }
    }
    
    # assure diagonal is all 1
    diag(sm) <- 1
    return(sm)  
}
