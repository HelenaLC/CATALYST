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

adaptCompensationSpillmat <- function(input_sm, output_names){
  y=input_sm
  check_spillMat(y)
  n = length(output_names)
  nms = output_names
  ms <- gsub("[[:alpha:][:punct:]]", "", nms)
  ff_chs <- nms[!is.na(ms)]
  sm_chs <- rownames(y)
  sm_cols <- colnames(y)
  y <- make_symetric(y)
  
  # check which channels of input flowFrame are not 
  # contained in spillover matrix and give warning

  
  # add them into the matrix
  sm <- matrix(diag(n), n, n, dimnames=list(nms, nms))
  sl_sm_cols <- sm_cols[sm_cols %in% ff_chs]
  sm[sl_sm_cols, sl_sm_cols] <- y[sl_sm_cols, sl_sm_cols]
  
  new_chs <- ff_chs[!ff_chs %in% sm_chs]
  if (length(new_chs) != 0) {
    # TODO is this now also caputuring cases where new masses could receive spillover?
    old_and_new_ms <- prep_spillMat(new_chs, sm_chs)
    old_ms <- old_and_new_ms$old_ms
    new_ms <- old_and_new_ms$new_ms
  } else {
    old_ms <- as.numeric(gsub("[[:punct:][:alpha:]]", "", sm_chs))
    new_ms <- new_ms
  }
  
  test <- (length(new_chs) != 0) && (any(inds <- old_ms %in% new_ms))
  if (test) {
    # check if any new masses were already present in the old masses
    # and add them to receive spillover according to the old masses
    
    # get the channels that correspond to the old_masses 
    # that have an aditional metal with the same weight
    y_col <- sm_chs[inds]
    names(y_col) <- as.character(old_ms[inds])
    # get all columns that are part of the affected masses
    fil <- ms %in% old_ms[inds]
    sm_col <- nms[fil]
    sm_col_ms <- as.character(ms[fil])
    # add the spillover
    sm[rownames(y), sm_col] <- y[, y_col[sm_col_ms]]
    for (m in unique(sm_col_ms)){
      mfil <- ms == m
      # set the spillover between channels of the same mass to 0
      # otherwise the linear system can get singular
      # diagonal elements will be set to 1 again later on
      sm[mfil, mfil] <- 0
    }
  }
  
  # check which channels of spillover matrix are missing in flowFrame
  # and drop corresponding rows and columns
  ex <- rownames(sm)[!rownames(sm) %in% ff_chs]
  if (length(ex) != 0)
    sm <- sm[!rownames(sm) %in% ex, !colnames(sm) %in% ex]
  
  # assure diagonal is all 1
  diag(sm) <- 1
  return(sm)  
}
