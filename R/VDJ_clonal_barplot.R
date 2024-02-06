#' Function to create stacked barplots to visualize clonal expansion per repertoire directly from a VDJ matrix (either from the minimal_VDJ() or VDJ_GEX_matrix())
#' @description Function to create stacked barplots to visualize clonal expansion per repertoire directly from a VDJ matrix (either from the minimal_VDJ() or VDJ_GEX_matrix()).
#' @param VDJ VDJ matrix (either from the minimal_VDJ() or VDJ_GEX_matrix()
#' @param counts.to.use The column name in the VDJ matrix of the clonotypes you want to use. Defaults to "clonotype_id".
#' @param group.by The column name in the VDJ matrix on which you want to seperate the repertoire plots. If the entire VDJ matrix is one repertoire, this argument should be "none" or empty.
#' @param expanded.colors Character vector. Colors to use for expanded clones. Should be more than 3 for better visibility. Defaults to a "darkorchid3"-based palette.
#' @param non.expanded.color Character. Color to use for non expanded clones. Defaults to "black"
#' @return Returns a list with a ggplot for each group.by element.
#' @export
#' @examples
#' out <- VDJ_expansion_stackedbarplot(VDJ, counts.to.use = "clonotype_id_10x",group.by = "sample_id")
#' do.call("grid.arrange", c(out, ncol = 6))
VDJ_clonal_barplot <- function(VDJ,
                               counts.to.use,
                               group.by,
                               expanded.colors,
                               non.expanded.color){
  
  if(missing(counts.to.use)) counts.to.use <- "clonotype_id"
  if(missing(group.by)) group.by <- "none"
  if(!(group.by %in% colnames(VDJ)) & group.by != "none"){stop("group.by is not in the VDJ object. Please provide a group.by that is in the column names of the VDJ object.")}
  if(missing(expanded.colors)) expanded.colors <- c("darkorchid4","darkorchid1","mediumorchid1","mediumpurple3")
  if(missing(non.expanded.color)) non.expanded.color <- "black"
  
  output_list <- list()
  
  #When all cells should be plotted in the same bar (group.by = "none") add a ubiquitous identifier
  if(group.by == "none"){
    VDJ$none <- "all_cells"
  }
  
  #Create existing combinations between counts.to.use and group.by
  VDJ$clonotype_group <- paste0(VDJ[,counts.to.use],"_",VDJ[,group.by])
  #Create a dataframe with the counts.to.use and group.by combinations and their frequency
  clonal_df <- as.data.frame(table(VDJ[,"clonotype_group"]),stringsAsFactors = FALSE)
  #Split the group.by from the counts.to.use
  clonal_df[,c("Var1","sample")] <- stringr::str_split_fixed(as.character(clonal_df$Var1), "_", 2)
  
  #Create a plot per group
  unordered_df <- clonal_df
  for(group in unique(unordered_df$sample)){
    clonal_df <- unordered_df[unordered_df$sample == group,]
    #Order the dataframe from most to least expanded
    clonal_df <- dplyr::arrange(clonal_df, desc(Freq))
    #Assign a color to each clonotype. Black when it is unexpanded, repeating expanded.colors for expanding clones.
    clonal_df[clonal_df$Freq == 1, "Var1"] <- "unexpanded"
    clonal_df$manual_colors <- rep_len(expanded.colors, length.out = nrow(clonal_df))
    clonal_df[clonal_df$Freq == 1,"manual_colors"] <- non.expanded.color
    #remove the first NA row
    clonal_df <- na.omit(clonal_df)
    #create a named vector to store the colors
    cl <- stats::setNames(clonal_df$manual_colors, clonal_df$Var1)
    #Create stacked barplot grouped on "sample" and colored on the frequency
    p <- ggplot2::ggplot(clonal_df, aes(x=as.factor(sample), y=Freq, fill=factor(Var1, levels = unique(Var1)))) +
      geom_bar(position = "fill", stat = "identity") +
      scale_fill_manual(values = cl)  +
      theme_minimal() +
      theme(text = element_text(size = 20),
            axis.text.y = element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()) +
      guides(fill = FALSE)
    #Add the ggplot to the output list
    output_list[[group]] <- p
  }
  #return a list where each element is a ggplot object
  return(output_list)
}


