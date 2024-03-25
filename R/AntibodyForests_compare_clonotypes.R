#' Function to compare tree topology of B cell lineages
#' @description Function to compare trees of clonotypes.
#' @param input - list - An AntibodyForests-object
#' @param min.nodes - integer - The minimum number of nodes in a tree to include in the comparison
#' @param distance.method - string - The method to calculate distance (default ...)
#' 'none'           : No distance metric, analyze similarity directly from distance.metrics
#' 'euclidean'      : 
#' 'jensen-shannon' : Jensen-Shannon distance between spectral density profiles of trees.
#' @param distance.metrics - string - If distance.method is "none" or "euclidean", these metrics will be used to calculate clusters and PCA/MDS dimensions and are used for plotting. (Default is mean.depth and nr.nodes)
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
#' @param clustering.method - string - Method to cluster trees (default none)
#' 'none'           : No clustering
#' 'mediods'        : Clustering based on the k-mediods method. The number of clusters is estimated based on the optimum average silhouette.
#' @param visualization.methods - string - The methods to analyze similarity (default PCA)
#' 'PCA'            : Scatterplot of the first two principal components. This is usefull when distance.method is "none".
#' 'MDS'            : Scatterplot of the first two dimensions using multidimensional scaling. Usefull for all distance methods
#' 'heatmap'        : A (clustered) heatmap of the distance between clonotypes. If distance.method is "none", euclidean distance will be calculated.
# #' @param metrics.to.visualize - string - Other metrics from the input to use for visualization.
#' @param plot.label - boolean - Label clonotypes in the PCA/MDS plot (default FALSE)
#' @param parallel If TRUE, the metric calculations are parallelized (default FALSE)
#' @param num.cores Number of cores to be used when parallel = TRUE (Defaults to all available cores - 1)
#' @return - list - Returns a distance matrix, clustering, and various plots based on visualization.methods
#' @export

AntibodyForests_compare_clonotypes <- function(input,
                                    min.nodes,
                                    distance.method,
                                    distance.metrics,
                                    clustering.method,
                                    visualization.methods,
                                    #metrics.to.visualize,
                                    plot.label,
                                    parallel,
                                    num.cores){
  
  #1. Set defaults and check for missing or incorrect input
  if(missing(input)){stop("Please provide a valid input object.")}
  if (class(input) != "AntibodyForests"){stop("The input is not in the correct format.")}
  if(missing(visualization.methods)){visualization.methods = "PCA"}
  #if(missing(metrics.to.visualize)){metrics.to.visualize = "none"}
  if(missing(distance.metrics)){distance.metrics = c("mean.depth", "nr.nodes")}
  if(missing(clustering.method)){clustering.method = "none"}
  if(missing(plot.label)){plot.label = F}
  if(missing(parallel)){parallel <- F}
  if(parallel == TRUE && missing(num.cores)){num.cores <- parallel::detectCores() -1}

  #3. Define functions
  calculate_euclidean <- function(df){
    #save names
    names <- rownames(df)
    #make all columns numeric
    df <- apply(df, 2, as.numeric)
    #Add rownames
    rownames(df) <- names
    #calculate euclidean distance on the scaled dataframe
    distance_matrix <- stats::dist(scale(df), method = "euclidean", diag = T, upper = T)

    return(distance_matrix)
  }
  
  #Calculate principle components
  calculate_PC <- function(df, to.scale){
    #Make sure df is a matrix
    df <- as.matrix(df)
    #Save names for later
    names <- rownames(df)
    #make all columns numeric
    df <- apply(df, 2, as.numeric)
    #Run a PCA and save the PCs in a dataframe
    pca_results <- as.data.frame(stats::prcomp(df, scale. = to.scale)$x)
    #Keep the first two PCs
    pca_results <- pca_results[,1:2]
    colnames(pca_results) <- c("Dim1", "Dim2")
    #Add tree names to the dataframe
    pca_results$tree <- names
    return(pca_results)
  }
  
  #Multidimensional scaling
  calculate_MDS <- function(df, distance.method){
    names <- rownames(df)
    #If distance method is "none", a (euclidean) distance structure should first be computed
    if(distance.method == "none"){
      distance_matrix <- calculate_euclidean(df)
    }else{distance_matrix <- df}
    #Compute classical metric multidimensional scaling
    results <- as.data.frame(stats::cmdscale(distance_matrix))
    colnames(results) <- c("Dim1", "Dim2")
    #Add tree names to the dataframe
    results$tree <- names
    
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
    distance_matrix <- RPANDA::JSDtree(phylo = phylo_list, meth = c("normal1"))
    return(distance_matrix)
  }
  
  cluster_mediods <- function(df, distance.method){
    #If distance method is "none", a (euclidean) distance structure should first be computed
    if(distance.method == "none"){
      distance_matrix <- calculate_euclidean(distance_matrix)
    }else{distance_matrix <- df}
    
    #Define the max number of clusters
    max_cluster <- dim(as.matrix(distance_matrix))[1]-1
    #Perform clustering
    mediods <- fpc::pamk(distance_matrix,krange=1:max_cluster)
    clusters <- mediods$pamobject$clustering
    #Assign the clusters to the tree names
    names(clusters) <- labels(as.matrix(distance_matrix))[[1]]
    return(clusters)
  }
  
  plot <- function(df, color, name, plot.label){
    p <- ggplot2::ggplot(df, ggplot2::aes(x=Dim1,y=Dim2, color=.data[[color]])) +
      ggplot2::geom_point(size=5) +
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(size = 20)) +
      ggplot2::ggtitle(name)
    
    if(plot.label){
      p <- p + ggrepel::geom_label_repel(ggplot2::aes(label = tree))
    }
    
    return(p)
  }
  
  #Plot a heatmap of the distance matrix
  heatmap <- function(df, clusters){
    #If there are clusters make clustered heatmap
    if(!(is.null(clusters))){
      #Make df of the clusters
      clusters <- as.data.frame(clusters)
      colnames(clusters) <- "cluster"
      clusters$cluster <- as.factor(clusters$cluster)
      #Plot the clustered heatmap
      p <- pheatmap::pheatmap(as.matrix(df),
                              annotation_col = clusters)
    }
    #If there are no clusters
    else{p <- pheatmap::pheatmap(as.matrix(df))}
    
    return(p)
  }

  #Create empty list to store the output
  output <- list()
  
  #1. Calculate metrics
  metric_df <- AntibodyForests_metrics(input = input, min.nodes = min.nodes, metrics = distance.metrics,
                          parallel = parallel, num.cores = num.cores)
  
  #1. Distance matrix
  #If distance.method is "none", metric dataframe will be used directly for visualization
  if (distance.method == "none"){
    #Row with NA will be removed for visualization
    distance_matrix <- stats::na.omit(metric_df)
    }
  #If distance.method is "jensen-shannon" use the trees with at least min.nodes in the AntibodyForests-object to calculate distance
  else if (distance.method == "jensen-shannon"){
    distance_matrix <- calculate_JS(input, min.nodes)
  }
  else if(distance.method == "euclidean"){
    #Row with NA will be removed for visualization
    metric_df <- stats::na.omit(metric_df)
    distance_matrix <- calculate_euclidean(metric_df)
  }
  #Store distance matrix (or metric dataframe when distance.method is "none") in the final output
  output[["distance_matrix"]] <- distance_matrix
  
  
  #2. Dimension reduction
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
  
  #3. Clustering
  if (clustering.method != "none"){
    if (clustering.method == "mediods"){
      clusters <- cluster_mediods(distance_matrix, distance.method)
    }
    
    output[["clustering"]] <- clusters
  }
  
  #4. Plotting
  #Create empty list to store plots
  plot_list <- list()
  for (method in visualization.methods){
    #Plot the default and store in plot list
    if(method == "heatmap" && clustering.method == "none"){
      plot_list[[paste0(method, "_default")]] <- heatmap(distance_matrix, cluster = NULL)
    }else if(method != "heatmap"){
      plot_list[[paste0(method, "_default")]] <- plot(results_list[[method]], color = "tree", name = method, plot.label) 
    }
    
    #Plot the clustering
    if (clustering.method != "none"){
      if(method == "heatmap"){
        plot_list[[paste0(method, "_clusters")]] <- heatmap(distance_matrix, cluster = clusters)
      }else{
        #Connect the clustering to the dimensions
        cluster_df <- cbind(results_list[[method]], cluster = as.factor(output[["clustering"]]))
        #Plot the clusters and store in plot list
        plot_list[[paste0(method, "_clusters")]] <- plot(cluster_df, color = "cluster", name = method, plot.label)
        
        #Exchange spectral.density for specific density metrics
        if ("spectral.density" %in% distance.metrics){
          #remove spectral.density
          distance.metrics <- distance.metrics[distance.metrics != "spectral.density"]
          #Add density specific metrics
          distance.metrics <- c(distance.metrics, "spectral.peakedness","spectral.asymmetry","spectral.principal.eigenvalue","modalities")
        }
        #Plot colored on optional metrics
        for (metric in distance.metrics){
          #Add this metric to the output dataframe
          metric_df <- stats::na.omit(metric_df)
          results_list[[method]][metric] <- as.numeric(metric_df[,metric])
          #Add extra plot to the output list
          plot_list[[paste0(method, "_", metric)]] <- plot(results_list[[method]], color = metric,
                                                           name = method, plot.label)
        }
      }
    }
  }
  
  #Add the plot list to the output
  output[["plots"]] <- plot_list
  
  #return output
  return(output)

}