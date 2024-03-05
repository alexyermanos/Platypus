#' Function to compare trees.
#' @description Function to compare trees.
#' @param input If within.clonotypes is FALSE: a list with an AntibodyForests-object (input[[1]]) and a metrics-dataframe (input[[2]]). If within.clonotypes is TRUE: a list with a list of AntibodyForests-objects (input[[1]]) and a list of metrics-dataframes (input[[2]])
#' @param within.clonotypes If TRUE: input should contain multiple metric dataframes (default FALSE)
#' @param min.nodes The minimum number of nodes in a tree to include in the comparison
#' @param distance.method The method to calculate distance (default ...)
#' 'none'           : No distance metric, analyze similarity directly from tree metrics in the input dataframe(s)
#' 'euclidean'      : Euclidean distance of the number of edges for each leaf (depth) between different trees of the same clonotype.
#' 'jensen-shannon' : Jensen-Shannon distance between spectral density profiles of trees.
#' @param distance.metrics If distance.method is "none", these metrics from the metric dataframe will be used to calculate PCA/MDS dimensions.
#' @param visualization.methods The methods to analyze similarity (default PCA)
#' 'PCA'            : Scatterplot of the first two principal components. This is usefull when distance.method is "none".
#' 'MDS'            : Scatterplot of the first two dimensions using multidimensional scaling. Usefull for all distance methods
#' @param metrics.to.visualize Other metrics from the input to use for visualization.
#' @return If within.clonotypes = TRUE: Returns a list of distance matrices for each clonotype and various plots based on visualization.methods and metrics.to.visualize. If within.clonotypes = FALSE returns a single distance matrix and various plots based on visualization.methods and metrics.to.visualize.
#' @export

AntibodyForests_compare <- function(input,
                                    within.clonotypes,
                                    min.nodes,
                                    distance.method,
                                    distance.metrics,
                                    visualization.methods,
                                    metrics.to.visualize){
  
  #1. Set defaults and check for missing input
  if(missing(input)){stop("Please provide a valid input object.")}
  if(missing(within.clonotypes)){within.clonotypes = F}
  if(missing(visualization.methods)){visualization.methods = "PCA"}
  if(missing(metrics.to.visualize)){metrics.to.visualize = "none"}
  if(missing(distance.metrics)){distance.metrics = "all"}
  
  #2. Check if the input is correct
  if ((length(input) != 2) ||
      (within.clonotypes == F && class(input[[1]][[1]][[1]][["igraph"]]) != "igraph" && 
       !("matrix" %in% class(input[[2]]))) ||
      (within.clonotypes == T && class(input[[1]][[1]][[1]][[1]][["igraph"]]) != "igraph" &&
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
  if((all(metrics.to.visualize != "none")) &&
     ((within.clonotypes == F && !(all(metrics.to.visualize %in% colnames(input[[2]])))) ||
      (within.clonotypes == T && !(all(metrics.to.visualize %in% lapply(input[[2]], colnames)))))){
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
    #make al columns numeric
    df <- apply(df, 2, as.numeric)
    #Run a PCA and save the PCs in a dataframe
    pca_results <- as.data.frame(stats::prcomp(df, scale. = to.scale)$x)
    #Keep the first two PCs
    pca_results <- pca_results[,1:2]
    colnames(pca_results) <- c("Dim1", "Dim2")
    #Add tree names to the dataframe
    pca_results$tree <- rownames(pca_results)
    
    return(pca_results)
  }
  
  #Multidimensional scaling
  calculate_MDS <- function(df, distance.method){
    #If distance method is "none", a (euclidean) distance structure should first be computed
    if(distance.method == "none"){
      distance_matrix <- stats::dist(df, method = "euclidean", diag = T, upper = T)
    }else{distance_matrix <- df}
    #Compute classical metric multidimensional scaling
    results <- as.data.frame(stats::cmdscale(distance_matrix))
    colnames(results) <- c("Dim1", "Dim2")
    #Add tree names to the dataframe
    results$tree <- rownames(results)
    
    return(results)
  }
  
  
  calculate_JS <- function(af, min.nodes){
    #Convert the igraph trees to phylo trees and store in list
    phylo_list <- list()
    for(sample in names(af)){
      for(clonotype in names(af[[sample]])){
        #Only keep trees with a minimum number of nodes (min.nodes)
        if (igraph::vcount(af[[sample]][[clonotype]][['igraph']]) >= min.nodes){
          phylo_tree <- AntibodyForests_phylo(af[[sample]][[clonotype]][['igraph']], solve_multichotomies = F)
          phylo_list[[paste0(sample,".",clonotype)]] <- phylo_tree
        }
      }
    }
    #Calculate distance 
    distance_matrix <- RPANDA::JSDtree(phylo = phylo_list, meth = c("normal2"))
    return(distance_matrix)
  }
  
  plot <- function(df, color, name){
    p <- ggplot2::ggplot(df, ggplot2::aes(x=Dim1,y=Dim2, color=.data[[color]])) +
      ggplot2::geom_point(size=5) +
      # ggrepel::geom_label_repel(aes(label = tree))+
      # ggrepel::geom_label_repel(data = . %>% dplyr::mutate(label = ifelse(tree %in% c("reclonotyped.clonotype41", "reclonotyped.clonotype46"),
      #                                                              tree, "")),
      #                           aes(label = label))+
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
      if (all((distance.metrics != "all"))){distance_matrix <- distance_matrix[,distance.metrics]}
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
    #Create empty list to store dimension reduction results
    results_list <- list()
    
    if ("PCA" %in% visualization.methods){
      #Calculate the principal components
      results <- calculate_PC(distance_matrix, to.scale = T)
      #Store results
      results_list[["PCA"]] <- results
      }
    
    if ("MDS" %in% visualization.methods){
      #Calculate the principal components
      results <- calculate_MDS(distance_matrix, distance.method)
      #Store results
      results_list[["MDS"]] <- results
      }
    
    for (method in visualization.methods){
      #Plot the default and store in plot list
      plot_list[[paste0(method, "_default")]] <- plot(results_list[[method]], color = "tree", name = method)
      
      #Plot colored on optional metrics
      if (all(metrics.to.visualize != "none")){
        for (metric in metrics.to.visualize){
          #Add this metric to the output dataframe
          metric_df <- stats::na.omit(input[[2]])
          results_list[[method]][metric] <- as.numeric(metric_df[,metric])
          #Add extra plot to the output list
          plot_list[[paste0(method, "_", metric)]] <- plot(results_list[[method]], color = metric,
                                         name = method)
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
        
        if (all(metrics.to.visualize != "none")){
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