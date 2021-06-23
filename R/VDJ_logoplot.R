#' Plots a logoplot of the CDR3 aminoacid region ! For platypus v3 and more flexibility use VDJ_logoplot_vector
#' @param VDJ.object Output of the VDJ_analyze function
#' @param length_cdr3 Integer indicating the length of the CDR3 regions that are selected to be plotted.
#' @return Returns the logo plot.
#' @export
#' @examples
#' \dontrun{
#' VDJ_logoplot <- function(VDJ.object = VDJ_analyze.out,length_cdr3 = 10)
#' }
VDJ_logoplot <- function(VDJ.object,
                         length_cdr3) {
  require(ggseqlogo)
  if (missing(length_cdr3) | length_cdr3 < min(stringr::str_length(VDJ.object[[1]]$CDR3_aa_pasted))) {
    length_cdr3 <- min(stringr::str_length(VDJ.object[[1]]$CDR3_aa_pasted))
  }
  if (length_cdr3 > max(stringr::str_length(VDJ.object[[1]]$CDR3_aa_pasted))) {
    length_cdr3 <- max(stringr::str_length(VDJ.object[[1]]$CDR3_aa_pasted))
  }
  print(ggseqlogo::ggseqlogo(VDJ.object[[1]]$CDR3_aa_pasted[which(nchar(VDJ.object[[1]]$CDR3_aa_pasted)==length_cdr3)], method='prob', seq_type="aa"))
}
