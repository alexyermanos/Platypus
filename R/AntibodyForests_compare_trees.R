#' Function to compare different trees from the same clonotype using the euclidean distance between the germline-to-node depth of all nodes.
#' @description Function to compare different trees from the same clonotype. For example to compare various graph construction and phylogenetic reconstruction methods.
#' @param input A list of AntibodyForests-objects as output from the function AntibodyForests(). These objects should contain the same samples/clonotypes. Please name the objects in the list according to their tree-construction method.
#' @param min.nodes The minimum number of nodes in a tree to include in the comparison, this includes the germline. Default is 2 (this includes all trees).
#' @param depth - string - Method to calculate the germline-to-node depth (default edge.count)
#' 'edge.count'   : The number of edges between each node and the germline
#' 'edge.length'  : The sum of edge lengths between each node and the germline
#' @param clustering.method Method to cluster trees (default NULL)
#' NULL             : No clustering
#' 'mediods'        : Clustering based on the k-mediods method. The number of clusters is estimated based on the optimum average silhouette.
#' @param visualization.methods The methods to analyze similarity (default NULL)
#' NULL             : No visualization
#' 'PCA'            : Scatterplot of the first two principal components.
#' 'MDS'            : Scatterplot of the first two dimensions using multidimensional scaling.
#' "heatmap'        : Heatmap of the distance
#' @param parallel If TRUE, the depth calculations are parallelized across clonotypes (default FALSE)
#' @param num.cores Number of cores to be used when parallel = TRUE. (Defaults to all available cores - 1)
#' @return A list with all clonotypes that pass the min.nodes threshold including the distance matrix, possible clustering and visualization
#' @export

AntibodyForests_compare_trees <- function(input,
                                    min.nodes,
                                    depth,
                                    clustering.method,
                                    visualization.methods,
                                    parallel,
                                    num.cores){
  
  #1. Set defaults and check for missing input
  if(missing(input)){stop("Please provide a valid input object.")}
  if(missing(min.nodes)){min.nodes = 2}
  if(missing(depth)){depth = "edge.count"}
  if(missing(visualization.methods)){visualization.methods = NULL}
  if(missing(clustering.method)){clustering.method = NULL}
  if(missing(parallel)){parallel <- FALSE}
  # If 'parallel' is set to TRUE but 'num.cores' is not specified, the number of cores is set to all available cores - 1
  if(parallel == TRUE && missing(num.cores)){num.cores <- parallel::detectCores() -1}
  
  #2. Check if the input is correct
  if (class(input[[1]][[1]][[1]][["igraph"]]) != "igraph"){stop("The input is not in the correct format.")}
  if(!(all(equal <- sapply(input, function(x) all.equal(names(x), names(input[[1]])))))){stop("The input list does not contain the same clonotypes.")}
  if(min.nodes > max(unlist(lapply(input, function(x){lapply(x,function(y){lapply(y, function(y){lapply(y$nodes, function(z){z$size})})})})))){
    stop("min.nodes is larger than the biggest clonotype.")}
  if(!(depth %in% c("edge.count", "edge.length"))){stop("Unvalid depth.")}
  if(!(is.null(clustering.method)) && clustering.method != "mediods"){stop("Unvalid clustering method.")}
  if(!(is.null(visualization.methods)) && all(!(visualization.methods %in% c("PCA", "MDS", "heatmap")))){stop("Unvalid visualization methods.")}
  
  
  #3. Define functions
  #Calculate the number of edges between certain nodes and the germline of a single tree
  calculate_depth <- function(tree, nodes, depth){
    if (depth == "edge.count"){
      #Get the shortest paths between each node and the germline
      paths <- igraph::shortest_paths(tree, from = "germline", to = nodes, output = "both")
      #Set names to the list of vpath
      names(paths$epath) <- names(unlist(lapply(paths$vpath, function(x){tail(x,n=1)})))
      #Reorder according to node number
      paths$epath <- paths$epath[paste0("node",sort(as.numeric(stringr::str_sub(names(paths$epath), start = 5))))]
      #Get the number of edges along the shortes paths
      depths <- unlist(lapply(paths$epath, length))
      #Set names to the vector of edge counts
      names(depths) <- names(paths$epath)
    }
    if (depth == "edge.length"){
      #Get the total length of shortest paths between each node and the germline
      depths <- igraph::distances(tree, v = "germline", to = nodes, algorithm = "dijkstra",
                                    weights = as.numeric(igraph::edge_attr(tree)$edge.length))
      #Reorder according to node number
      depths <- depths[,paste0("node",sort(as.numeric(stringr::str_sub(colnames(depths), start = 5))))]
    }
    #Return the named vector of depths per node
    return(depths)
  }
  
  #Remove entries from (nested) lists that contain NA
  remove_na_entries <- function(lst) {
    if (!is.list(lst)) {
      return(lst)
    } else {
      cleaned <- lapply(lst, remove_na_entries)
      cleaned <- cleaned[sapply(cleaned, function(x) !all(is.na(x)))]
      return(cleaned)
    }
  }
  
  #Calculate the euclidean distance between node depth for each clonotype in the depth_list
  calculate_euclidean <- function(depth_list){
    distance_list <- list()
    samples <- unique(unlist(lapply(depth_list, names)))
    for (sample in samples){
      #Get the clonotype names in this sample
      clonotypes <- names(depth_list[[1]][[sample]])
      for (clonotype in clonotypes){
        #Create a dataframe where the rows are the tree construction methods and the columns are the depth per node
        nodes <- names(depth_list[[1]][[sample]][[clonotype]])
        depth_df <- matrix(ncol = length(nodes), nrow = 0)
        colnames(depth_df) <- nodes
        for (method in names(depth_list)){depth_df <- rbind(depth_df, depth_list[[method]][[sample]][[clonotype]])}
        rownames(depth_df) <- names(depth_list)
        
        #Calculate the euclidean distance between the tree methods for this clonotype
        euclidean_matrix <- stats::dist(depth_df, method = "euclidean", diag = T, upper = T)
        
        #Add distance matrix to the list
        distance_list[[paste0(sample,".",clonotype)]] <- euclidean_matrix
      }
    }
    return(distance_list)
  }
  
  cluster_mediods <- function(df){
    distance_matrix <- df
    #Define the max number of clusters
    max_cluster <- dim(as.matrix(distance_matrix))[1]-1
    #Perform clustering
    mediods <- fpc::pamk(distance_matrix,krange=1:max_cluster)
    clusters <- mediods$pamobject$clustering
    #Assign the clusters to the tree names
    names(clusters) <- labels(df)
    return(clusters)
  }
  
  #Calculate principle components
  calculate_PC <- function(df, to.scale){
    names <- labels(df)
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
  calculate_MDS <- function(df){
    names <- labels(df)
    #Compute classical metric multidimensional scaling
    results <- as.data.frame(stats::cmdscale(df))
    colnames(results) <- c("Dim1", "Dim2")
    #Add tree names to the dataframe
    results$tree <- names
    
    return(results)
  }
  
  #Plot the first two dimensions of a PCA or MDS
  plot <- function(df, color, name){
    p <- ggplot2::ggplot(df, ggplot2::aes(x=Dim1,y=Dim2, color=as.factor(.data[[color]]))) +
      ggplot2::geom_point(size=5) +
      ggrepel::geom_label_repel(ggplot2::aes(label = tree))+
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(size = 20)) +
      ggplot2::ggtitle(name)
    
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
  

  #4. Calculate the depths
  # If 'parallel' is set to TRUE, the depth calculation is parallelized across the trees
  if(parallel){
    # Retrieve the operating system
    operating_system <- Sys.info()[['sysname']]
    # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
    if(operating_system %in% c('Linux', 'Darwin')){
        #Go over each tree in clonotype and create a depth list
        depth_list <- parallel::mclapply(mc.cores = num.cores, input, function(object){
          parallel::mclapply(mc.cores = num.cores, object, function(sample){
            parallel::mclapply(mc.cores = num.cores, sample, function(clonotype){
              #Only calculate depth for trees with at least min.nodes
              if (igraph::vcount(clonotype$igraph) >= min.nodes){
                #Calculate the depth for all nodes except the germline
                depth_vector <- calculate_depth(clonotype$igraph,
                                                nodes = igraph::V(clonotype$igraph)[names(igraph::V(clonotype$igraph)) != "germline"],
                                                depth)
                return(depth_vector)
              }else{
                return(NA)
              }
              
            })
          })
        })
    }
    # If the operating system is Windows, "parLapply" is used for parallelization
    if(operating_system == "Windows"){
      # Create cluster
      cluster <- parallel::makeCluster(num.cores)
      #Go over each tree in each of the AntibodyForests objects and create a list of depths
      depth_list <- parallel::parLapply(cluster, input, function(object){
        parallel::parLapply(cluster, object, function(sample){
          parallel::parLapply(cluster, sample, function(clonotype){
            #Only calculate depth for trees with at least min.nodes
            if (igraph::vcount(clonotype$igraph) >= min.nodes){
              #Calculate the depth for all nodes except the germline
              depth_vector <- calculate_depth(clonotype$igraph,
                                              nodes = igraph::V(clonotype$igraph)[names(igraph::V(clonotype$igraph)) != "germline"],
                                              depth)
              return(depth_vector)
            }else{
              return(NA)
            }
          })
        })
      })
      # Stop cluster
      parallel::stopCluster(cluster)
    }
  }
  # If 'parallel' is set to FALSE, the network inference is not parallelized
  if(!parallel){
    #Go over each tree in each of the AntibodyForests objects and create a list of depths
    depth_list <- lapply(input, function(object){
      lapply(object, function(sample){
        lapply(sample, function(clonotype){
          #Only calculate depth for trees with at least min.nodes
          if (igraph::vcount(clonotype$igraph) >= min.nodes){
            #Calculate the depth for all nodes except the germline
            depth_vector <- calculate_depth(clonotype$igraph,
                                            nodes = igraph::V(clonotype$igraph)[names(igraph::V(clonotype$igraph)) != "germline"],
                                            depth)
            return(depth_vector)
          }else{
            return(NA)
          }
        })
      })
    })
  }
  #Remove NA (clonotypes with less then min.nodes nodes)
  depth_list <- remove_na_entries(depth_list)
  
  #5. Calculate the euclidean distance between trees for each clonotype per sample
  output_list <- calculate_euclidean(depth_list)

  #7. For each clonotype
  output_list <- lapply(output_list, function(distance_matrix){
    #Create inner list
    temp_list <- list()
    #Add the distance matrix to this list
    temp_list[["distance.matrix"]] <- distance_matrix
    
    #Clustering
    if(!(is.null(clustering.method))){
      #K-mediods clustering
      if (clustering.method == "mediods"){
        #Get clusters
        clusters <- cluster_mediods(distance_matrix)
      }
      #Add to the list
      temp_list[["clusters"]] <- clusters
    }
    
    #Visualization
    if(!(is.null(visualization.methods))){
      #PCA analysis
      if("PCA" %in% visualization.methods){
        #Get PCA dimensions
        pca <- calculate_PC(distance_matrix, to.scale = F)
        #If there are clusters calculated
        if(!(is.null(clustering.method))){
          #Add clusters to the pca dataframe
          pca <- cbind(pca, "clusters" = clusters)
          #Create plot
          plot <- plot(pca, color = "clusters", name = "PCA")
        }
        #If there are no clusters
        else{
          #Create plot
          plot <- plot(pca, color = "tree", name = "PCA")
        }
        #Add to the list
        temp_list[["PCA"]] = plot
      }
      #MDS analysis
      if("MDS" %in% visualization.methods){
        #Only do MDS when there are more than 2 trees to compare
        if (ncol(as.matrix(distance_matrix)) > 2){
          #Get MDS dimensions
          mds <- calculate_MDS(distance_matrix)
          #If there are clusters calculated
          if(!(is.null(clustering.method))){
            #Add clusters to the pca dataframe
            mds <- cbind(mds, "clusters" = clusters)
            #Create plot
            plot <- plot(mds, color = "clusters", name = "MDS")
          }
          #If there are no clusters
          else{
            #Create plot
            plot <- plot(mds, color = "tree", name = "MDS")
          }
          #Add to the list
          temp_list[["MDS"]] = plot
        }else{
          temp_list[["MDS"]] = "Need at least 3 trees to compute MDS"
        }
      }
      #Heatmap
      if("heatmap" %in% visualization.methods){
        #If there are clusters calculated
        if(!(is.null(clustering.method))){
          #Create plot
          hm <- heatmap(distance_matrix, clusters = clusters)
        }
        #If there are no clusters
        else{
          #Create plot
          hm <- heatmap(distance_matrix, clusters = NULL)
        }
        #Add to the list
        temp_list[["Heatmap"]] <- hm
      }
      
    }
    
    return(temp_list)
  })
  
  return(output_list)


}