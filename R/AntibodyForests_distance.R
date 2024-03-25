#' Function to make a grouped boxplot of distance between nodes from specific groups and the germline of lineage trees constructed with AntibodyForests.
#' @description Function to compare trees.
#' @param input AntibodyForests-object with node features to compare distance to the germline
#' @param distance - string - How to calculate the distance to the germline.
#' 'node.depth'     : Average of the sum of nodes on the shortest parth between germline and nodes from this group. (Default)
#' 'edge.length'    : Average of the sum of edge length of the shortest path between germline and nodes from this group.
#' @param min.nodes The minimum number of nodes for a tree to be included in this analysis (this included the germline)
#' @param groups Which groups to compare. These groups need to be in the node features of the AntibodyForests-object.
#' If you want to compare IgM and IgG for example, groups should be c("IgM, "IgG") (not "Isotypes")
#' @param colors Optionally specific colors for the group (Will be matched to the groups/names on alphabetical order).
#' @param text.size Font size in the plot (default 20).
#' @param parallel If TRUE, the metric calculations are parallelized across clonotypes. (default FALSE)
#' @export

AntibodyForests_distance <- function(input,
                                     distance,
                                     min.nodes,
                                     groups, 
                                     colors,
                                     text.size,
                                     parallel){
  
  #Set defaults and check for missing input
  if(missing(input)){stop("Please provide an AntibodyForests-object as input.")}
  if(missing(groups)){stop("Please provide groups to compare.")}
  if(missing(colors)){colors = scales::hue_pal()(length(groups))}
  if(missing(distance)){distance = "node.depth"}
  if(missing(text.size)){text.size = 20}
  if(missing(min.nodes)){min.nodes = 0}
  if(missing(parallel)){parallel <- F}
  #Check if group are in the metric dataframe
  #if(!(all(groups %in% colnames(metric_df)))){stop("Groups are not in the column names of the metric dataframe.")}
  
  #Calculate the average distance to the germline per group
  metric_df <- AntibodyForests_metrics(input,
                                       parallel = parallel,
                                       min.nodes = min.nodes,
                                       metrics = paste0("group.",distance))
  
  #Error if zero or only one tree is in the metric_df
  if(is.null(nrow(metric_df))){stop("Your AntibodyForests-object does not have enough trees that pass the min.nodes threshold.")}
  
  #Get the column of the groups to compare
  df <- as.data.frame(metric_df[,paste0(groups,".node.depth")])
  
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
    ggplot2::geom_line(ggplot2::aes(group=clonotype)) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = text.size),
                   legend.position = "none")  +
    ggplot2::scale_x_discrete(breaks=paste0(groups,".node.depth"),
                     labels=groups) +
    ggplot2::ggtitle(paste0("Distance (", distance, ") to germline"))

}
