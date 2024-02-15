#' Function to compare trees.
#' @description Function to compare trees.
#' @param input Output from AntibodyForests_metrics(), either one dataframe or a list of dataframes
#' @param within.clonotypes If TRUE: input should contain multiple metric dataframes (default FALSE)
#' @param distance.methods The metrics to be calculated (default ...)
#' euclidean'    : Euclidean distance of the number of edges for each leaf (depth) between different trees of the same clonotype.
#' @param similarity.methods The methods to analyze similarity (default PCA)
#' @param metrics.to.visualize Other metrics from the input to use for visualization.
#' @return Returns a list of distance matrices for each clonotype. Optionally various plots based on similarity.methods and metrics.to.visualize
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
  if(class(input) == "list" && length(unique(lapply(input, rownames))) != 1){
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
  if(metrics.to.visualize != "none" && !(metrics.to.visualize %in% lapply(input, colnames)[[1]])){
    stop("metrics.to.visualize need to be in the input metrics dataframes.")
  }
  
  #Functions to calculate distance
  calculate_euclidean <- function(tree_list, clonotype){
    #Create a matrix where each row is a tree and each column is the depth per node
    depth_matrix <- matrix(data = NA, ncol = as.numeric(tree_list[[1]][clonotype,"nodes"]))
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
  
  if ("PCA" %in% similarity.methods){
    #Save names of the trees
    list_names <- names(distance_list)
    
    #Create empty list to store plots
    plot_list <- list()
    
    #Do a PCA and plotting for each clonotype
    distance_list <- lapply(seq_along(distance_list), function(clonotype){
      #Calculate principle components
      pca_results <- as.data.frame(prcomp(distance_list[[clonotype]])$x)
      pca_results$tree <- rownames(pca_results)
      
      #Plot the PC1 and PC2 of the distance between the trees
      pca_plot_default <- ggplot2::ggplot(pca_results, aes(x=PC1,y=PC2,color=tree, label = tree)) +
        ggplot2::geom_point(size=5) +
        ggrepel::geom_label_repel()+
        ggplot2::theme_minimal() +
        ggplot2::theme(text = element_text(size = 20),
              legend.position = "none") +
        ggplot2::ggtitle(names(distance_list)[[clonotype]])
      plot_list[["default"]] <- pca_plot_default
      
      if (metrics.to.visualize != "none"){
        for (metric in metrics.to.visualize){
          #Add this metric to the PCA output dataframe
          pca_results[metric] <- NA
          for (tree in names(within_metrics)){
            pca_results[tree,metric] <- as.numeric(within_metrics[[tree]][names(distance_list)[[clonotype]], metric])
          }
          
          #Plot the PC1 and PC2 of the distance between the trees colored on metric
          pca_plot_metric <- ggplot2::ggplot(pca_results, aes(x=PC1,y=PC2, color=.data[[metric]], label = tree)) +
            ggplot2::geom_point(size=5) +
            ggrepel::geom_label_repel() +
            ggplot2::theme_minimal() +
            ggplot2::theme(text = element_text(size = 20)) +
            ggplot2::ggtitle(names(distance_list)[[clonotype]])
          
          #Add extra PCA plot to the output list
          plot_list[[metric]] <- pca_plot_metric
        }
      }
      #Add the distance matrix and the default PCA plot to the output list
      output_list <- list(distance_list[[clonotype]], plot_list)
      
      return(output_list)
    })
    
    #Name the clonotypes in the list
    names(distance_list) <- list_names
    
  }
  
  return(distance_list)
}