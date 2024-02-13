#' Function to calculate metrics for each tree in an AntibodyForests-object
#' @description Function to calculate metrics for each tree in an AntibodyForests-object
#' @param input AntibodyForests-object(s), output from AntibodyForests()
#' @param multiple.objects If TRUE: input should contain multiple AntibodyForests-objects (default FALSE)
#' @param metrics The metrics to be calculated (default mean.depth)
#' 'mean.depth'       : Mean of the number of edges connecting each node to the germline
#' 'all.depth'        : Number of edges connecting each node to the germline
#' 'sackin.index'     : Sum of the number of nodes between each node and the germline
#' 'colless.number'   : Sum of the absolute difference between the number of left- and right descendants for each node (this requires a tree to be binary!)
#' @return Returns a dataframe where the rows are trees and the columns are metrics
#' @export

AntibodyForests_metrics <- function(input,
                                    multiple.objects,
                                    metrics){
  
  #Stop when no input is provided
  if(missing(input)){stop("Please provide a valid input object.")}
  #Set defaults
  if(missing(multiple.objects)){multiple.objects = F}
  if(missing(metrics)){metrics <- "mean.depth"}
  #Check if the input is in the correct format.
  #If multiple.objects is TRUE, multiple AntibodyForests-objects should be in the input list, where the third item in the nested AntibodyForest-object should be of class "igraph"
  if((multiple.objects == F && class(input[[1]][[1]][[3]]) != "igraph")|| (multiple.objects == T && class(input[[1]][[1]][[1]][[3]]) != "igraph")){
    stop("The input is not in the correct AntibodyForests-object format.")}
  
  #Functions to calculate metrics
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
  calculate_metrics <- function(clonotype, metrics){
    
    #Create empty vector to store metrics
    metrics_vector <- c()
    
    if ("mean.depth" %in% metrics){
      depth <- calculate_mean_depth(clonotype$lineage.tree)
      metrics_vector["mean_depth"] <- depth
    }
    
    if ("all.depth" %in% metrics){
      all_depth <- calculate_all_depth(clonotype$lineage.tree)
      metrics_vector["all_depth"] <- all_depth
    }
    
    if ("sackin.index" %in% metrics){
      si <- calculate_sackin_index(clonotype$lineage.tree)
      metrics_vector["sackin_index"] <- si
    }
    
    # Traditional Colless Number on bifurcating phylo tree
    # if ("colless.number" %in% metrics){
    #   colless <- calculate_colless_number(clonotype$phylo.tree)
    #   metrics_vector["colless_number"] <- colless
    # }
   
    return(metrics_vector)
  }
  
  #Create a list of metric dataframes if there are multiple AntibodyForests-objects in the input
  if(multiple.objects == T){
    #Go over each tree in each of the AntibodyForests objects and create a metric list
    metric_list <- lapply(input, function(object){
      t(as.data.frame(lapply(object, function(sample){
        lapply(sample, function(clonotype){
          calculate_metrics(clonotype, metrics)
        })
      })))
    })
    
    return(metric_list)
  }
  #Create a single metric dataframe if there is only one AntibodyForests-object in the input
  else if(multiple.objects == F){
    #Go over each tree in the AntibodyForests object and create a metric dataframe
    metric_df <- t(as.data.frame(lapply(input, function(sample){
      lapply(sample, function(clonotype){
        calculate_metrics(clonotype, metrics)
      })
    })))
    
    return(metric_df)
  }

  
  
}





