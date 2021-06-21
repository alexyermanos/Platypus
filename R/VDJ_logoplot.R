#' Plots a logoplot of the CDR3 aminoacid region
#' @param VDJ.object Output of the VDJ_analyze function
#' @param length_cdr3 Integer indicating the length of the CDR3 regions that are selected to be plotted.
#' @return Returns the logo plot.
#' @export
#' @examples
#' \dontrun{
#' GEX_fraction_cells_per_cluster(automate_GEX.output, per_sample))
#' }
VDJ_logoplot <- function(VDJ.object,
                         length_cdr3) {
  require(ggseqlogo)
  if (missing(length_cdr3) | length_cdr3 < min(str_length(VDJ.object[[1]]$CDR3_aa_pasted))) {
    length_cdr3 <- min(str_length(VDJ.object[[1]]$CDR3_aa_pasted))
  }
  if (length_cdr3 > max(str_length(VDJ.object[[1]]$CDR3_aa_pasted))) {
    length_cdr3 <- max(str_length(VDJ.object[[1]]$CDR3_aa_pasted))
  }
  print(ggseqlogo(VDJ.object[[1]]$CDR3_aa_pasted[which(nchar(VDJ.object[[1]]$CDR3_aa_pasted)==length_cdr3)], method='prob', seq_type="aa"))
}
