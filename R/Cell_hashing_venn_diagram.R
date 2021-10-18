#' Creates Venn diagrams that quantify the number of clones that contain cells with one or multiple hashing barcodes.
#' @param clonotype_phylo.dataframes Output of the VDJ_clonotype function. Needs to be in the phylo.dataframes form (see example below). Specify in the VDJ_clonotype function if output should be per sample or global.
#' @param hashing_BC  Character vector. Provide the names of the hashing bracodes used. 
#' @param text_size Numeric. Size of the labels showing the count for each hashing barcode. Default is sets to 3.
#' @param show_percentage Logical. Size of the labels showing the count for each hashing barcode. Default is TRUE.
#' @param digits Numeric. Number of digits after the decimal point for percentages. Default is set to 1.
#' @param stroke_linetype Numeric. Line width of the circle borders. Default is set to 0.
#' @param set_name_size Numeric. Size of the labels describing the barcodes. Default is set to 4.
#' @param fill_alpha Numeric. Values between 0 an 1. Modify color transparency. Default is set to 0.7.
#' @return returns a list containing Venn diagrams plots. 
#' @export
#' @examples
#' \dontrun{
#' clonotype_phylo.dataframes <- VDJ_clonotype(VDJ.GEX.matrix = VDJ.matrix,VDJ.VJ.1chain = T, global.clonotype = F,clone.strategy="cdr3.aa", output.format = "phylo.dataframes")
#' Cell_hashing_venn_diagram(clonotype_phylo.dataframes = clonotype_phylo.dataframes, hashing_BC = hashing_BC)

#'}


Cell_hashing_venn_diagram <- function(clonotype_phylo.dataframes,hashing_BC,text_size,show_percentage,digits,stroke_linetype,set_name_size,fill_alpha){
  
  #define default parameters if not specified 
  if(missing(hashing_BC)) stop("Provide valid hashing barcode names")
  
  
  if(missing(text_size)) text_size = 3
  if(missing(show_percentage)) show_percentage = TRUE
  if(missing(digits)) digits = 1
  if(missing(stroke_linetype)) stroke_linetype = 0
  if(missing(set_name_size)) set_name_size = 4
  if(missing(fill_alpha)) fill_alpha = 0.7
  
  outlist <- list()
  plot_clonotypes_per_BC <- list()
  clonotypes_per_BC <- list()
  for (k in 1:length(clonotype_phylo.dataframes)){
    clonotypes_per_BC[[k]] <- list()
  }
  
  for (k in 1:length(clonotype_phylo.dataframes)){
    for (j in 1:length(hashing_BC)){
      clonotypes_per_BC[[k]][[j]] <- vector()
    }
    for (i in 1:length(clonotype_phylo.dataframes[[k]])){
      clonotype_phylo.dataframes[[k]][[i]]$FB_assignment <- substr(clonotype_phylo.dataframes[[k]][[i]]$FB_assignment,4,nchar(clonotype_phylo.dataframes[[k]][[i]]$FB_assignment))
      clonotype_phylo.dataframes[[k]][[i]]$FB_assignment[is.na(clonotype_phylo.dataframes[[k]][[i]]$FB_assignment)] <- "Not assignable"
      clonotype_phylo.dataframes[[k]][[i]]$FB_assignment[clonotype_phylo.dataframes[[k]][[i]]$FB_assignment == " assignable"] <- "Not assignable"
      for (j in 1:length(hashing_BC)){
        if (hashing_BC[j] %in% unique(clonotype_phylo.dataframes[[k]][[i]]$FB_assignment)){
          clonotypes_per_BC[[k]][[j]] <- append(clonotypes_per_BC[[k]][[j]],names(clonotype_phylo.dataframes[[k]])[i])
        }
      }
    }
    names(clonotypes_per_BC[[k]]) <- hashing_BC
    plot_clonotypes_per_BC[[k]] <- ggvenn::ggvenn(clonotypes_per_BC[[k]],text_size = text_size,show_percentage = show_percentage,digits = digits,stroke_linetype = stroke_linetype,set_name_size = set_name_size,fill_alpha = fill_alpha,fill_color = scales::hue_pal()(length(hashing_BC))) + ggplot2::ggtitle(clonotype_phylo.dataframes[[k]][[1]]$sample_id[1]) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  
  outlist <- plot_clonotypes_per_BC
  return(plot_clonotypes_per_BC)
}#stop function






