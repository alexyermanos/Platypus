#' Function to plot the average number of nodes between nodes from specific groups and the germline of lineage trees constructed with AntibodyForests.
#' @description Function to compare trees.
#' @param metric_df Dataframe with trees metrics (output from AntibodyForests_metrics())
#' @param groups Which groups to compare. These groups need to columns in the metric_df (and supplied as node.features in AntibodyForests())
#' @param names Optionally strings to rename the groups. Order of names needs to be the same.
#' @param colors Optionally specific colors for the group (Will be matched to the groups/names on alphabetical order).
#' @param text.size Font size in the plot (default 20).
#' @export

AntibodyForests_distance <- function(metric_df, 
                                     groups, 
                                     names,
                                     colors,
                                     text.size){
  
  #Set defaults and check for missing input
  if(missing(metric_df)){stop("Please provide a metric dataframe as input.")}
  if(missing(groups)){stop("Please provide groups to compare.")}
  if(missing(names)){names = NA}
  if(missing(colors)){colors = scales::hue_pal()(length(groups))}
  if(missing(text.size)){text.size = 20}
  #Check if group are in the metric dataframe
  if(!(all(groups %in% colnames(metric_df)))){stop("Groups are not in the column names of the metric dataframe.")}
  
  #Get the column of the groups to compare
  df <- as.data.frame(metric_df[,groups])
  
  #Optionally names
  if(all(!(is.na(names)))){
    colnames(df) <- names
  }
  
  #Add clonotype as column
  df$clonotype <- rownames(df)
  
  #Only keep clonotypes that have nodes of all groups
  df <- as.data.frame(na.omit(df))
  #Check if there are clonotypes left after NA removal
  if(nrow(df) == 0){stop("No trees contain nodes from all groups.")}
  
  #Transform dataframe for visualization
  df <- tidyr::pivot_longer(df, cols=colnames(df)[1:ncol(df)-1],
                             names_to='group',
                             values_to='depth')
  
  #Plot the grouped boxplots with lines
  ggplot2::ggplot(df, ggplot2::aes(x=group, y=depth, fill=group)) + 
    ggplot2::geom_boxplot()+ 
    ggplot2::geom_point()+ 
    ggplot2::scale_fill_manual(values=colors) +
    ggplot2::geom_line(aes(group=clonotype)) +
    theme_classic() +
    ggplot2::theme(text = element_text(size = text.size))  +
    ggtitle("Node distance to germline")

}
