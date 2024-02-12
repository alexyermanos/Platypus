#' Function to calculate metrics for each tree in an AntibodyForests-object
#' @description Function to calculate metrics for each tree in an AntibodyForests-object
#' @param AF AntibodyForests-object, output from AntibodyForests()
#' @param metrics The metrics to be calculated
#' 'mean.depth'       : Mean of the number of edges connecting each node to the germline
#' 'all.depth'        : Number of edges connecting each node to the germline
#' 'sackin.index'     : Sum of the number of nodes between each node and the germline
#' 'colless.number'   : Sum of the absolute difference between the number of left- and right descendants for each node (this requires a tree to be binary!)
#' @return Returns a dataframe where the rows are trees and the columns are metrics
#' @export

AntibodyForests_metrics <- function(AF,
                                    metrics){
  
  calculate_mean_depth <- function(tree){
    paths <- igraph::shortest_paths(tree, from = "germline", output = "both")
    depth <- mean(unlist(lapply(paths$epath, length)))
    return(depth)
  }
  
  calculate_all_depth <- function(tree){
    paths <- igraph::shortest_paths(tree, from = "germline", output = "both")
    depths <- paste(unlist(lapply(paths$epath, length)), collapse = "_")
    return(depths)
  }
  
  calculate_sackin_index <- function(tree){
    paths <- igraph::shortest_paths(tree, from = "germline", output = "both")
    depth <- sum(unlist(lapply(paths$vpath, length)))
    return(depth)
  }
  
  calculate_colless_number <- function(tree){
    colless <- phyloTop::colless.phylo(ape::as.phylo(tree), normalise = T)
    return(colless)
  }
  
  #Calculate the metrics for a tree
  calculate_metrics <- function(tree, metrics){
    
    #Create empty vector to store metrics
    metrics_vector <- c()
    
    if ("mean.depth" %in% metrics){
      depth <- calculate_mean_depth(tree)
      metrics_vector["mean_depth"] <- depth
    }
    
    if ("all.depth" %in% metrics){
      all_depth <- calculate_all_depth(tree)
      metrics_vector["all_depth"] <- all_depth
    }
    
    if ("sackin.index" %in% metrics){
      si <- calculate_sackin_index(tree)
      metrics_vector["sackin_index"] <- si
    }
    
    if ("colless.number" %in% metrics){
      colless <- calculate_colless_number(tree)
      metrics_vector["colless_number"] <- colless
    }
   
    return(metrics_vector)
  }

  #Go over each tree in the AntibodyForests object and create a metric dataframe
  metric_df <- t(as.data.frame(lapply(AF, function(sample){
    lapply(sample, function(clonotype){
      calculate_metrics(clonotype$lineage.tree, metrics)
    })
  })))
  
  return(metric_df)
}





