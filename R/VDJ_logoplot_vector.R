#' Plots a logoplot of the CDR3 aminoacid region
#' @param cdr3.vector A character vector of aa sequences. This is to increase flexibility of this function. Such a sequence vector may be retrieved from the VDJ_analyse function output on a clonotype level or from the avdjgexders_assemble function output on a per cell level. Additionally, any length of sequence may be used (e.g. HCDR3 only or H and LCDR3 pasted together)
#' @param length_cdr3 Integer or character. Defaults to "auto". Sets the length of the CDR3 regions that are selected to be plotted. If set to auto, the most frequently appearing length in the vector will be used
#' @param seq_type passed to ggseqlogo. Can be set to "aa", "dna", "rna" or "other" 
#' @return Returns the logo plot.
#' @export
#' @examples
#' \dontrun{
#' GEX_fraction_cells_per_cluster(automate_GEX.output, per_sample))
#' }

VDJ_logoplot_vector <- function(cdr3.vector, 
                                length_cdr3,
                                seq_type) {
  require(ggseqlogo)
  if(missing(length_cdr3)) length_cdr3 <- "auto"
  if(missing(seq_type)) seq_type <- "auto"
  
  cdr3.vector <- as.character(cdr3.vector)
  cdr3.vector[is.na(cdr3.vector)] <- ""
  if("" %in% cdr3.vector){
  cdr3.vector <- cdr3.vector[-c(which(cdr3.vector == ""))]}

  if(length_cdr3 == "auto"){
    length_cdr3 <- as.numeric(names(which.max(table(nchar(cdr3.vector)))))
    print("Length distribution of input sequences")
    print(table(nchar(cdr3.vector)))
  } else if (length_cdr3 < min(nchar(cdr3.vector))) {
    length_cdr3 <- min(nchar(cdr3.vector))
  } else if (length_cdr3 > max(nchar(cdr3.vector))) {
    length_cdr3 <- max(nchar(cdr3.vector))
  }
  if(!as.character(length_cdr3) %in% names(table(nchar(cdr3.vector)))){
    print("No CDR3s of chosen length present. Setting length_cdr3 to mode of length of input sequences")
    length_cdr3 <- as.numeric(names(which.max(table(nchar(cdr3.vector)))))
  }
  print(paste0("Returning logoplot based on ", length(cdr3.vector[which(nchar(cdr3.vector)==length_cdr3)]), " sequences of a total ", length(cdr3.vector), " input sequences"))
  return(ggseqlogo(cdr3.vector[which(nchar(cdr3.vector)==length_cdr3)], method='prob', seq_type= seq_type))
}
