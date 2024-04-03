#' Function to make a grouped boxplot of metrics from clusters of clonotypes
#' @description Function to compare metrics between clusters of clontoypes
#' @param af - list - AntibodyForests-object with node features to compare distance to the germline
#' @param clusters - named integer - The clusters as output from AntibodyForests_compare_clonotypes
#' @param metrics - string - The metrics to be calculated per tree
#' 'nr.nodes'         : The total number of nodes
#' 'nr.cells'         : The total number of cells in this clonotype
#' 'mean.depth'       : Mean of the number of edges connecting each node to the germline
#' 'mean.edge.length' : Mean of the edge lengths between each node and the germline
#' 'group.depth'      : Mean of the number of edges connecting each node per group (node.features of the AntibodyForests-object) to the germline. (default FALSE)
#' 'sackin.index'     : Sum of the number of nodes between each node and the germline
#' 'spectral.density' : Metrics of the spectral density profiles (calculated with package RPANDA)
#'    - peakedness            : Tree balance
#'    - asymmetry             : Shallow or deep branching events
#'    - principal eigenvalue  : Phylogenetic diversity
#'    - modalities            : The number of different structures within the tree
#' @param min.nodes The minimum number of nodes for a tree to be included in this analysis (this included the germline). This should be the same as for the AntibodyForests_compare_clonotypes() functions.
#' @param colors - string -  Optionally specific colors for the clusters
#' @param text.size Font size in the plot (default 20).
#' @param parallel If TRUE, the metric calculations are parallelized across clonotypes. (default FALSE)
#' @return - list - A list with boxplots per metric
#' @export

AntibodyForests_clusters <- function(af,
                                     clusters,
                                     metrics,
                                     min.nodes,
                                     colors,
                                     text.size,
                                     parallel){
  
  #Set defaults and check for missing input
  if(missing(af)){stop("Please provide an AntibodyForests-object as input.")}
  if(missing(clusters)){stop("Please provide clusters as input.")}
  if(missing(metrics)){stop("Please provide metrics to calculate.")}
  if(missing(min.nodes)){min.nodes = 0}
  if(missing(colors)){colors = scales::hue_pal()(length(unique(clusters)))}
  if(missing(text.size)){text.size = 20}
  if(missing(parallel)){parallel <- F}
  #Check if group are in the metric dataframe
  #if(!(all(groups %in% colnames(metric_df)))){stop("Groups are not in the column names of the metric dataframe.")}
  
  #Calculate the metrics
  metric_df <- AntibodyForests_metrics(af,
                                       parallel = parallel,
                                       min.nodes = min.nodes,
                                       metrics = metrics)
  
  #Check if clonotypes are the same between metric_df and clusters
  if (nrow(metric_df) < length(clusters)){stop("Make sure that min.nodes threshold is not higher then min.nodes used when running AntibodyForests_compare_clonotypes().")}
  if (!(all(rownames(clusters) %in% names(metric_df)))){stop("The names of the clonotypes in the AntibodyForests-object and the clusters should be the same.")}
  
  #Add clusters to the metric_df
  cluster_df <- merge(metric_df, as.data.frame(clusters), by = "row.names")
  
  #Get the calculated metrics
  metrics <- colnames(cluster_df)[!(colnames(cluster_df) %in% c("Row.names", "clusters"))]
  
  #Make a plot for each metric and store in a list
  output_list <- list()
  for (metric in metrics){
    #Plot the grouped boxplots
    p <- ggplot2::ggplot(cluster_df, ggplot2::aes(x=as.factor(clusters), y=as.numeric(.data[[metric]]), fill=as.factor(clusters))) + 
      ggplot2::geom_boxplot() + 
      ggplot2::geom_jitter(color="black", size=1) +
      ggplot2::scale_fill_manual(values=colors) +
      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = text.size),
                     legend.position = "none") +
      ggplot2::xlab("Cluster") +
      ggplot2::ggtitle(paste0(metric, " per cluster"))
    
    #Add to output list
    output_list[[metric]] <- p
  }
  
  #return output list with plots
  return(output_list)
}