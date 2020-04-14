# ==============================================================================
# Estimate separation cutoffs
#' @rdname data
#' @name data
#' @aliases 
#' raw_data sample_ff sample_key ss_exp mp_cells isotope_list 
#' PBMC_fs PBMC_panel PBMC_md merging_table
#' @title Example data sets
#' 
#' @description 
#' \itemize{
#' \item{Concatenation & Normalization\describe{
#' \item{\code{raw_data}}{
#' a \code{\link{flowSet}} with 3 experiments, each containing 2'500 raw 
#' measurements with a variation of signal over time. Samples were mixed 
#' with DVS beads capture by mass channels 140, 151, 153, 165 and 175.}}}
#' 
#' \item{Debarcoding\describe{
#' \item{\code{sample_ff}}{
#' a \code{\link{flowFrame}} following a 6-choose-3 barcoding scheme 
#' where mass channels 102, 104, 105, 106, 108, and 110 were used for labeling 
#' such that each of the 20 individual barcodes are positive for exactly 3 
#' out of the 6 barcode channels.}
#' \item{\code{sample_key}}{
#' a \code{data.frame} of dimension 20 x 6 with sample names as row and 
#' barcode masses as column names. Contains a binary code of length 6 for each
#' sample in \code{sample_ff}, e.g. 111000, as its unique identifier.}}}
#' 
#' \item{Compensation\describe{
#' \item{\code{ss_exp}}{
#' a \code{\link{flowFrame}} with 20'000 events. 
#' Contains 36 single-antibody stained controls where beads were stained 
#' with antibodies captured by mass channels 139, 141 through 156, and 
#' 158 through 176, respectively, and pooled together.}
#' \item{\code{mp_cells}}{
#' a \code{\link{flowFrame}} with 5000 spill-affected
#' multiplexed cells and 39 measurement parameters.}
#' \item{\code{isotope_list}}{
#' a named list of isotopic compositions for all elements within 75 through 
#' 209 u corresponding to the CyTOF mass range at the time of writing.}}}
#' 
#' \item{Differential Analysis\describe{
#' \item{\code{PBMC_fs}}{ 
#' a \code{\link{flowSet}} with PBMCs samples from 6 patients. For each sample, 
#' the expression of 10 cell surface and 14 signaling markers was measured 
#' before (REF) and upon BCR/FcR-XL stimulation (BCRXL) with B cell receptor/
#' Fc receptor crosslinking for 30', resulting in a total of 12 samples.}
#' \item{\code{PBMC_panel}}{
#' a 2 column \code{data.frame} that contains each marker's
#' column name in the FCS file, and its targeted protein marker.}
#' \item{\code{PBMC_md}}{
#' a \code{data.frame} where each row corresponds to a sample, 
#' and with columns describing the experimental design.}
#' \item{\code{merging_table}}{
#' a 20 x 2 table with "old_cluster" IDs and "new_cluster" labels
#' to exemplify manual cluster merging and cluster annotation.}}}
#' }
#' 
#' @return 
#' see descriptions above.
#' 
#' @examples
#' ### example data for concatenation & normalization:
#'     # raw measurement data
#'     data(raw_data)
#'   
#' ### example data for debarcoding:
#'     # 20 barcoded samples
#'     data(sample_ff)
#'     # 6-choose-3 barcoding scheme
#'     data(sample_key)
#' 
#' ### example data for compensation:
#'     # single-stained control samples
#'     data(ss_exp)
#'     # multiplexed cells
#'     data(mp_cells)
#'     
#' ### example data for differential analysis:
#'     # REF vs. BCRXL samples
#'     data(PBMC_fs)
#'     # antigen panel & experimental design
#'     data(PBMC_panel, PBMC_md)
#'     # exemplary manual merging table
#'     data(merging_table)
#' 
#' @references
#' Bodenmiller, B., Zunder, E.R., Finck, R., et al. (2012). 
#' Multiplexed mass cytometry profiling of cellular states perturbed by 
#' small-molecule regulators. \emph{Nature Biotechnology} \bold{30}(9): 858-67.
#' 
#' Coursey, J.S., Schwab, D.J., Tsai, J.J., Dragoset, R.A. (2015). 
#' Atomic weights and isotopic compositions, 
#' (available at http://physics.nist.gov/Comp).
#' 
#' @author Helena L Crowell \email{helena.crowell@@uzh.ch}

NULL