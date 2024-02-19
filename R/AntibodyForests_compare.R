#' Function to compare trees.
#' @description Function to compare trees.
#' @param input Output from AntibodyForests_metrics(), either one dataframe or a list of dataframes
#' @param within.clonotypes If TRUE: input should contain multiple metric dataframes (default FALSE)
#' @param distance.methods The metrics to be calculated (default ...)
#' 'none'         : No distance metric, analyze similarity directly from tree metrics in the input dataframe(s)
#' 'euclidean'    : Euclidean distance of the number of edges for each leaf (depth) between different trees of the same clonotype.
#' @param similarity.methods The methods to analyze similarity (default PCA)
#' @param metrics.to.visualize Other metrics from the input to use for visualization.
#' @return If within.clonotypes = TRUE: Returns a list of distance matrices for each clonotype and various plots based on similarity.methods and metrics.to.visualize. If within.clonotypes = FALSE returns only plots.
#' @export

AntibodyForests_compare <- function(input,
                                    within.clonotypes,
                                    distance.methods,
                                    similarity.methods,
                                    metrics.to.visualize){
  
  #Stop when no input is provided
  if(missing(input)){stop("Please provide a valid input object.")}
  #Set defaults
  if(missing(within.clonotypes)){within.clonotypes = F}
  if(missing(similarity.methods)){similarity.methods = "PCA"}
  if(missing(metrics.to.visualize)){metrics.to.visualize = "none"}
  #Check if the input is a list when comparing within clonotypes.
  if(within.clonotypes == T && class(input) != "list"){
    stop("Please provide a list of metric dataframes when comparing within clonotypes.")
  }
  #When comparing within clonotypes, check if the same trees are in the metric dataframes.
  if(within.clonotypes == T && class(input) == "list" && length(unique(lapply(input, rownames))) != 1){
    stop("Row names of the metric dataframes are not matching. Make sure to select the same 
         trees (in the same order) when comparing within clonotypes.")
  }
  #If euclidean distance is selected, check if comparing within clonotypes
  if(within.clonotypes == F && "euclidean" %in% distance.methods){
    stop("Euclidean distance can only be calculated between different trees from the same clonotype.")
  }
  #If euclidean distance is selected, check if all_depth is in the metrics dataframes
  if("euclidean" %in% distance.methods && !("all_depth" %in% lapply(input, colnames)[[1]])){
    stop("all_depth needs to be calculated in order to compare euclidean distance.")
  }
  #Check if all metrics in metrics.to.visualize are in the input dataframes
  if((within.clonotypes == T && metrics.to.visualize != "none" && !(metrics.to.visualize %in% lapply(input, colnames)[[1]])) ||
     (within.clonotypes == F && metrics.to.visualize != "none" && !(metrics.to.visualize %in% colnames(input)))){
    stop("metrics.to.visualize need to be in the input metrics dataframes.")
  }
  
  #Functions to calculate distance
  calculate_euclidean <- function(tree_list, clonotype){
    #Create a matrix where each row is a tree and each column is the depth per node
    depth_matrix <- matrix(data = NA, ncol = as.numeric(tree_list[[1]][clonotype,"nr_nodes"]))
    for(index in 1:length(tree_list)){
      all_depth <- tree_list[[index]][clonotype,"all_depth"]
      #split the all_depth string into a vector
      all_depth <- stringr::str_split(all_depth, "_")[[1]]
      depth_matrix <- rbind(depth_matrix, as.numeric(all_depth))
    }
    depth_matrix <- na.omit(depth_matrix)
    rownames(depth_matrix) <- names(tree_list)
    
    #Calculate the distance matrix from the depth matrix
    euclidean_matrix <- stats::dist(depth_matrix, method = "euclidean", diag = T, upper = T)
    
    return(euclidean_matrix)
  }
  
  #Calculate principle components
  calculate_PC <- function(df, to.scale){
    #Run a PCA and save the PCs in a dataframe
    pca_results <- as.data.frame(stats::prcomp(df)$x, scale. = to.scale)
    
    #Add tree names to the dataframe
    pca_results$tree <- rownames(pca_results)
    
    return(pca_results)
  }
  
  plot_PC <- function(df, color, name){
    p <- ggplot2::ggplot(df, aes(x=PC1,y=PC2, color=.data[[color]], label = tree)) +
      ggplot2::geom_point(size=5) +
      ggrepel::geom_label_repel()+
      ggplot2::theme_minimal() +
      ggplot2::theme(text = element_text(size = 20)) +
      ggplot2::ggtitle(name)
    
    return(p)
  }
  
  #1. Calculate distance matrix
  
  if (within.clonotypes == T){
    #Create an output list to store the distance matrix for each clonotype
    distance_list <- as.list(unique(lapply(input, rownames))[[1]])
    names(distance_list) <- unique(lapply(input, rownames))[[1]]
    
    #Calculate the distance for each clonotype
    for (clonotype in unique(lapply(input, rownames))[[1]]){
      if ("euclidean" %in% distance.methods){
        distance_matrix <- calculate_euclidean(input, clonotype)
        distance_list[[clonotype]] <- distance_matrix
      }
    }
  }
  
  #2. Similarity analysis and visualization
  
  #PCA on the distance matrix per clonotype when within.clonotypes is TRUE
  if ("PCA" %in% similarity.methods && within.clonotypes == T){
    #Save names of the trees
    list_names <- names(distance_list)
    
    #Do a PCA and plotting for each clonotype
    distance_list <- lapply(seq_along(distance_list), function(clonotype){
      #Create empty list to store plots
      plot_list <- list()
      
      #Calculate principle components
      pca_results <- calculate_PC(distance_list[[clonotype]], to.scale = F)
      
      #Plot the PC1 and PC2 of the distance between the trees and color on tree
      plot_list[["default"]] <- plot_PC(df = pca_results, color = "tree",
                                          name = names(distance_list)[[clonotype]])
      
      if (metrics.to.visualize != "none"){
        for (metric in metrics.to.visualize){
          #Add this metric to the PCA output dataframe
          pca_results[metric] <- NA
          for (tree in names(input)){
            pca_results[tree,metric] <- as.numeric(input[[tree]][names(distance_list)[[clonotype]], metric])
          }
          
          #Add extra PCA plot to the output list
          #Plot the PC1 and PC2 of the distance between the trees and color on metric
          plot_list[[metric]] <- plot_PC(psa_results, color = metric,
                                            name = names(distance_list)[[clonotype]])
        }
      }
      #Add the distance matrix and the default PCA plot to the output list
      output_list <- list(distance_list[[clonotype]], plot_list)
      
      return(output_list)
    })
    
    #Name the clonotypes in the list
    names(distance_list) <- list_names
    
    return(distance_list)
    
  }
  #PCA directly on the input on all clonotypes when within.clonotypes is FALSE
  if ("PCA" %in% similarity.methods && within.clonotypes == F){
    #Calculate principle components
    pca_results <- calculate_PC(input, to.scale = T)
    
    #Create empty list to store plots
    plot_list <- list()
    
    #Plot the PC1 and PC2 of the distance between the trees and color on tree
    plot_list[["default"]] <- plot_PC(pca_results, color = "tree", name = "All trees")
    
    if (metrics.to.visualize != "none"){
      for (metric in metrics.to.visualize){
        #Add this metric to the PCA output dataframe
        pca_results[metric] <- as.numeric(input[,metric])
        
        #Add extra PCA plot to the output list
        #Plot the PC1 and PC2 of the distance between the trees and color on metric
        plot_list[[metric]] <- plot_PC(pca_results, color = metric,
                                       name = "All trees")
    
      }
    }
    
    return(plot_list)
  }
  
  
}