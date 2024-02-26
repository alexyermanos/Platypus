#' Function to compare trees.
#' @description Function to compare trees.
#' @param input If within.clonotypes is FALSE: a list with an AntibodyForests-object (input[[1]]) and a metrics-dataframe (input[[2]]). If within.clonotypes is TRUE: a list with a list of AntibodyForests-objects (input[[1]]) and a list of metrics-dataframes (input[[2]])
#' @param within.clonotypes If TRUE: input should contain multiple metric dataframes (default FALSE)
#' @param min.nodes The minimum number of nodes in a tree to include in the comparison
#' @param distance.method The method to calculate distance (default ...)
#' 'none'           : No distance metric, analyze similarity directly from tree metrics in the input dataframe(s)
#' 'euclidean'      : Euclidean distance of the number of edges for each leaf (depth) between different trees of the same clonotype.
#' 'jensen-shannon' : Jensen-Shannon distance between spectral density profilesof trees.
#' @param visualization.methods The methods to analyze similarity (default PCA)
#' @param metrics.to.visualize Other metrics from the input to use for visualization.
#' @return If within.clonotypes = TRUE: Returns a list of distance matrices for each clonotype and various plots based on visualization.methods and metrics.to.visualize. If within.clonotypes = FALSE returns only plots.
#' @export

AntibodyForests_compare <- function(input,
                                    within.clonotypes,
                                    min.nodes,
                                    distance.method,
                                    visualization.methods,
                                    metrics.to.visualize){
  
  #1. Set defaults and check for missing input
  if(missing(input)){stop("Please provide a valid input object.")}
  if(missing(within.clonotypes)){within.clonotypes = F}
  if(missing(visualization.methods)){visualization.methods = "PCA"}
  if(missing(metrics.to.visualize)){metrics.to.visualize = "none"}
  
  #2. Check if the input is correct
  if ((length(input) != 2) ||
      (within.clonotypes == F && class(input[[1]][[1]][[1]][[3]]) != "igraph" && 
       !("matrix" %in% class(input[[2]]))) ||
      (within.clonotypes == T && class(input[[1]][[1]][[1]][[1]][[3]]) != "igraph" &&
       !("matrix" %in% class(input[[2]][[1]])))){
    stop("The input is not in the correct format.")}
  if(within.clonotypes == T && !(all(names(input[[2]]) %in% names(input[[1]])))){
    stop("The names of AntibodyForests-objects and metric-dataframes need to be the same.")
  }
  if((within.clonotypes == F && 
      !(all(unlist(lapply(strsplit(x = rownames(input[[2]]), split = "\\."), function(x){x[1]})) %in% names(input[[1]]))) && 
      !(all(unlist(lapply(strsplit(x = rownames(input[[2]]), split = "\\."), function(x){x[2]})) %in% lapply(input, names)))) ||
     (within.clonotypes == T && 
      !(all(unlist(lapply(strsplit(rownames(input[[2]][[1]]), "\\."), function(x){x[1]})) %in% names(input[[1]][[1]]))) && 
      !(all(unlist(lapply(strsplit(rownames(input[[2]][[1]]), "\\."), function(x){x[2]})) %in% lapply(input[[1]], names))))){
    stop("The trees in the metric-dataframe(s) need to be present in the AntibodyForests-object(s).")
  }
  if(within.clonotypes == F && distance.method == "euclidean"){
    stop("Euclidean distance can only be calculated between different trees from the same clonotype.")
  }
  if(distance.method == "euclidean" && !(all(lapply(input[[2]], function(x){"all_depth" %in% colnames(x)})))){
    stop("all_depth needs to be calculated in order to compare euclidean distance.")
  }
  if((metrics.to.visualize != "none") &&
     ((within.clonotypes == F && !(metrics.to.visualize %in% colnames(input[[2]]))) ||
      (within.clonotypes == T && !(metrics.to.visualize %in% lapply(input[[2]], colnames))))){
    stop("metrics.to.visualize need to be in the input metrics dataframe(s).")
  }

  #3. Define functions
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
  
  calculate_JS <- function(af, min.nodes){
    #Convert the igraph trees to phylo trees and store in list
    phylo_list <- list()
    for(sample in names(af)){
      for(clonotype in names(af[[sample]])){
        #Only keep trees with a minimum number of nodes (min.nodes)
        if (igraph::vcount(af[[sample]][[clonotype]][['lineage.tree']]) >= min.nodes){
          phylo_tree <- AntibodyForests_phylo(af[[sample]][[clonotype]][['lineage.tree']], solve_multichotomies = F)
          phylo_list[[paste0(sample,".",clonotype)]] <- phylo_tree
        }
      }
    }
    #Calculate distance 
    distance_matrix <- RPANDA::JSDtree(phylo = phylo_list, meth = c("standard"))
    return(distance_matrix)
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

  
  
  ##############
  if (within.clonotypes == F){
    #Create empty list to store the output
    output <- list()
    
    #1. Distance matrix
    #If distance.method is "none", metric dataframe will be used directly for visualization
    if (distance.method == "none"){
      #Row with NA will be removed for visualization
      distance_matrix <- stats::na.omit(input[[2]])
    }
    #If distance.method is "jensen-shannon" use the trees with at least min.nodes in the AntibodyForests-object to calculate distance
    else if (distance.method == "jensen-shannon"){
      distance_matrix <- calculate_JS(input[[1]], min.nodes)
    }
    #Store distance matrix (or metric dataframe when distance.method is "none") in the final output
    output[["distance_matrix"]] <- distance_matrix
    
    #2. Visualization
    #Create empty list to store plots
    plot_list <- list()
    
    if ("PCA" %in% visualization.methods){
      #Calculate the principal components
      pca_results <- calculate_PC(distance_matrix, to.scale = T)
      
      #Plot the PCA, color on tree (default) and store in plot list
      plot_list[["default"]] <- plot_PC(pca_results, color = "tree", name = "all")
      
      #Plot the PCA colored on optional metrics
      if (metrics.to.visualize != "none"){
        for (metric in metrics.to.visualize){
          #Add this metric to the PCA output dataframe
          metric_df <- stats::na.omit(input[[2]])
          pca_results[metric] <- as.numeric(metric_df[,metric])
          
          #Add extra PCA plot to the output list
          #Plot the PC1 and PC2 of the distance between the trees and color on metric
          plot_list[[metric]] <- plot_PC(pca_results, color = metric,
                                         name = "All trees")
          }
      }
    }
    
    #Add the plot list to the output
    output[["plots"]] <- plot_list
  }else if(within.clonotypes == T){
    #Create empty list to store the output
    output <- list()
    
    #1. Distance matrix
    #If distance.method is "none", metric dataframe will be used directly for visualization
    if (distance.method == "euclidean"){
      #Create an output list to store the distance matrix for each clonotype
      distance_list <- as.list(unique(lapply(input[[2]], rownames))[[1]])
      names(distance_list) <- unique(lapply(input[[2]], rownames))[[1]]
      
      #Calculate the distance for each clonotype
      for (clonotype in names(distance_list)){
          distance_matrix <- calculate_euclidean(input[[2]], clonotype)
          distance_list[[clonotype]] <- distance_matrix
        }
      
    }
    #Store distance matrices (or metric dataframe when distance.method is "none") in the final output
    output[["distance_matrix"]] <- distance_list
    
    #2. Visualization
    
    if ("PCA" %in% visualization.methods){
      #Save names of the trees
      list_names <- names(distance_list)
      
      #Do a PCA and plotting for each clonotype
      plot_list <- lapply(seq_along(distance_list), function(clonotype){
        #Create empty list to store plots
        temp_list <- list()
        
        #Calculate principle components
        pca_results <- calculate_PC(distance_list[[clonotype]], to.scale = F)
        
        #Plot the PC1 and PC2 of the distance between the trees and color on tree
        temp_list[["default"]] <- plot_PC(df = pca_results, color = "tree",
                                          name = names(distance_list)[[clonotype]])
        
        if (metrics.to.visualize != "none"){
          for (metric in metrics.to.visualize){
            #Add this metric to the PCA output dataframe
            pca_results[metric] <- NA
            for (tree in names(input[[2]])){
              metric_df <- stats::na.omit(input[[2]])
              pca_results[tree,metric] <- as.numeric(input[[2]][[tree]][names(distance_list)[[clonotype]], metric])
            }
            
            #Add extra PCA plot to the output list
            #Plot the PC1 and PC2 of the distance between the trees and color on metric
            temp_list[[metric]] <- plot_PC(psa_results, color = metric,
                                           name = names(distance_list)[[clonotype]])
          }
        }

        return(temp_list)})
      
      #Name the clonotypes in the list
      names(plot_list) <- list_names
    }
    #Add the plot list to the output
    output[["plots"]] <- plot_list 
  }  

  
  #return output
  return(output)
}