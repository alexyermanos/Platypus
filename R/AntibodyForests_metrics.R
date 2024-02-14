#' Function to calculate metrics for each tree in an AntibodyForests-object
#' @description Function to calculate metrics for each tree in an AntibodyForests-object
#' @param input AntibodyForests-object(s), output from AntibodyForests()
#' @param multiple.objects If TRUE: input should contain multiple AntibodyForests-objects (default FALSE)
#' @param metrics The metrics to be calculated (default mean.depth)
#' 'mean.depth'       : Mean of the number of edges connecting each node to the germline
#' 'all.depth'        : Number of edges connecting each node to the germline
#' 'sackin.index'     : Sum of the number of nodes between each node and the germline
#' 'colless.number'   : Sum of the absolute difference between the number of left- and right descendants for each node (this requires a tree to be binary!)
#' @param parallel If TRUE, the metric calculations are parallelized across AntibodyForests-objects. (default FALSE)
#' @param num.cores Number of cores to be used when parallel = TRUE. (Defaults to all available cores - 1)
#' @return Returns a dataframe where the rows are trees and the columns are metrics
#' @export

AntibodyForests_metrics <- function(input,
                                    multiple.objects,
                                    metrics,
                                    parallel,
                                    num.cores){
  
  #Stop when no input is provided
  if(missing(input)){stop("Please provide a valid input object.")}
  #Set defaults
  if(missing(multiple.objects)){multiple.objects = F}
  if(missing(metrics)){metrics <- "mean.depth"}
  if(missing(parallel)){parallel <- FALSE}
  #Check if the input is in the correct format.
  #If multiple.objects is TRUE, multiple AntibodyForests-objects should be in the input list, where the third item in the nested AntibodyForest-object should be of class "igraph"
  if((multiple.objects == F && class(input[[1]][[1]][[3]]) != "igraph")|| (multiple.objects == T && class(input[[1]][[1]][[1]][[3]]) != "igraph")){
    stop("The input is not in the correct AntibodyForests-object format.")}
  # If 'parallel' is set to TRUE but 'num.cores' is not specified, the number of cores is set to all available cores - 1
  if(parallel == TRUE && missing(num.cores)){num.cores <- parallel::detectCores() -1}
  
  #Functions to calculate metrics
  calculate_mean_depth <- function(tree){
    paths <- igraph::shortest_paths(tree, from = "germline", output = "both")
    depth <- mean(unlist(lapply(paths$epath, length)))
    return(depth)
  }
  
  calculate_all_depth <- function(tree){
    paths <- igraph::shortest_paths(tree, from = "germline", output = "both")
    #Set names to the list of vpath
    names(paths$epath) <- names(unlist(lapply(paths$vpath, function(x){tail(x,n=1)})))
    #Reorder the vpaths according to node number
    paths$epath <- paths$epath[c("germline",sort(names(paths$epath)[names(paths$epath) != "germline"]))]
    #Concatenate the depths in a string separated by "_"
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
  
  #Calculate the metrics for a clonotype
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
  


  # If 'parallel' is set to TRUE, the metric calculation is parallelized across the clonotypes
  if(parallel){
    
    # Retrieve the operating system
    operating_system <- Sys.info()[['sysname']]
    
    # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
    if(operating_system %in% c('Linux', 'Darwin')){
      #Create a list of metric dataframes if there are multiple AntibodyForests-objects in the input
      if(multiple.objects == T){
        #Go over each tree in each of the AntibodyForests objects and create a metric list
        metric_list <- parallel::mclapply(mc.cores = num.cores, input, function(object){
          t(as.data.frame(parallel::mclapply(mc.cores = num.cores, object, function(sample){
            parallel::mclapply(mc.cores = num.cores, sample, function(clonotype){
              calculate_metrics(clonotype, metrics)
            })
          })))
        })
        
        return(metric_list)
      }
      #Create a single metric dataframe if there is only one AntibodyForests-object in the input
      else if(multiple.objects == F){
        #Go over each tree in the AntibodyForests object and create a metric dataframe
        metric_df <- t(as.data.frame(parallel::mclapply(mc.cores = num.cores,input, function(sample){
          parallel::mclapply(mc.cores = num.cores,sample, function(clonotype){
            calculate_metrics(clonotype, metrics)
          })
        })))
        
        return(metric_df)
      }
    }
    
    # If the operating system is Windows, "parLapply" is used for parallelization
    if(operating_system == "Windows"){
      # Create cluster
      cluster <- parallel::makeCluster(num.cores)
      
      if(multiple.objects == T){
        #Go over each tree in each of the AntibodyForests objects and create a metric list
        metric_list <- parallel::parLapply(cluster, input, function(object){
          t(as.data.frame(parallel::parLapply(cluster, object, function(sample){
            parallel::parLapply(cluster, sample, function(clonotype){
              calculate_metrics(clonotype, metrics)
            })
          })))
        })
        # Stop cluster
        parallel::stopCluster(cluster)
        
        return(metric_list)
      }
      #Create a single metric dataframe if there is only one AntibodyForests-object in the input
      else if(multiple.objects == F){
        #Go over each tree in the AntibodyForests object and create a metric dataframe
        metric_df <- t(as.data.frame(parallel::parLapply(cluster,input, function(sample){
          parallel::parLapply(cluster,sample, function(clonotype){
            calculate_metrics(clonotype, metrics)
          })
        })))
        # Stop cluster
        parallel::stopCluster(cluster)
        
        return(metric_df)
      }

    }
    
  }
  
  # If 'parallel' is set to FALSE, the network inference is not parallelized
  if(!parallel){
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
  
  
}





