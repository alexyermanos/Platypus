#' Infer intermediate nodes in the minimum spanning trees/ sequences similiarity networks created by the AntibodyForests function

#'@description Intermediate nodes are expanded/inferred based on the edge weight between two existing nodes: for example, of node 1 and node 2 are connected by an edge of weight = 3, 2 new nodes are added in-between and all resulting edges have weight = 1.
#' @param trees AntibodyForests object/list of AntibodyForests objects - the resulting sequence similarity or minimum spanning tree networks from the AntibodyForests function.
#' @param parallel boolean - whether to execute the main subroutine in parallel or not. Requires the 'parallel' R package.


#' @return nested list of AntibodyForests objects or single AntibodyForests object, with the resulting networks having expanded/inferred intermediate nodes.
#' @examples
#' \dontrun{
#' AntibodyForests_expand_intermediates(trees)
#'}



AntibodyForests_expand_intermediates <- function(trees,
                                                 parallel
                                                ){

  if(missing(trees)) stop('Please input a nested list of AntibodyForests objects which you want to expand.')
  if(missing(parallel)) parallel <- T


  #SUBROUTINE 1: expands the intermediate nodes for a single tree/similarity network.
  expand_intermediates <- function(tree){

     adjacency_matrix <- as.matrix(igraph::as_adjacency_matrix(tree@tree, attr = 'weight'))
     network_df <- igraph::as_data_frame(tree@tree, what = 'vertices')
     directed <- igraph::is_directed(tree@tree)

     if(!directed){
       node_number <- sum(adjacency_matrix)/2 + nrow(adjacency_matrix) - length(which(adjacency_matrix!=0))/2
     }else{
       node_number <- sum(adjacency_matrix) + nrow(adjacency_matrix) - length(which(adjacency_matrix!=0))
     }

     added_intermediates <- node_number - nrow(adjacency_matrix)

     if(added_intermediates == 0){
       return(tree@tree)
     }

     final_adjacency_matrix <- matrix(0, node_number, node_number)
     final_adjacency_matrix[1:nrow(adjacency_matrix), 1:nrow(adjacency_matrix)] <- adjacency_matrix

     new_node <- nrow(adjacency_matrix) + 1
     for(i in 1:nrow(adjacency_matrix)){
       neighbours_to_expand <- which(adjacency_matrix[i, ]>1)
       if(length(neighbours_to_expand)==0){
         next
       }

       for(neighbour in neighbours_to_expand){
         edge_weight <- adjacency_matrix[i, neighbour]
         prev_added_node <- i

         while(edge_weight>1){
           final_adjacency_matrix[prev_added_node, new_node] <- 1
           if(!directed){
             final_adjacency_matrix[new_node, prev_added_node] <- 1
           }
           prev_added_node <- new_node
           new_node <- new_node + 1
           edge_weight <- edge_weight - 1
         }

        final_adjacency_matrix[i, neighbour] <- 0
        final_adjacency_matrix[neighbour, i] <- 0
        adjacency_matrix[i, neighbour] <- 0
        adjacency_matrix[neighbour, i] <- 0
        final_adjacency_matrix[new_node-1, neighbour] <- 1

        if(!directed){
          final_adjacency_matrix[neighbour, new_node-1] <- 1
        }
       }
     }

     clonotype_id <- paste0(unique(network_df$clonotype_id), collapse = ';')
     sample_id <- paste0(unique(network_df$sample_id), collapse = ';')


     #Add the intermediate rows into the network_df
     for(i in 1:added_intermediates){

       network_df <- rbind(network_df, NA)
       network_df$node_type[nrow(network_df)] <- 'intermediate'
       network_df$clonotype_id[nrow(network_df)] <- clonotype_id
       network_df$sample_id[nrow(network_df)] <- sample_id
       network_df$label[nrow(network_df)] <- nrow(network_df)
       network_df$cell_number[nrow(network_df)] <- 1

     }

     features <- tree@feature_names

     if(length(features)!=0){
       for(feature in features){
         network_df[feature][network_df$node_type == 'intermediate',] <- 'intermediate'
         network_df[paste0(feature,'_counts')][network_df$node_type == 'intermediate',] <- 1
       }
    }

     edgelist  <-  igraph::graph_from_adjacency_matrix(final_adjacency_matrix) %>%
                   igraph::as_data_frame()

     network_df$name <- 1:nrow(network_df)

     g <- igraph::graph_from_data_frame(edgelist, directed = directed, vertices = network_df)
     igraph::V(g)$label <- 1:length(igraph::V(g))
     igraph::V(g)$name <- 1:length(igraph::V(g))
     igraph::E(g)$weight <- 1
     return(g)
  }

  #MAIN LOOP
  if(inherits(trees, 'list')){
    for(i in 1:length(trees)){
      expanded_list <- vector(mode = "list", length = length(trees[[i]]))

      if(parallel){
        requireNamespace('parallel')
        cores <- parallel::detectCores()
        expanded_list <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                 FUN = function(x) {x %>% expand_intermediates()

                                                                    })

      }else{
        expanded_list <- lapply(trees[[i]], function(x) {x %>% expand_intermediates()
                                                              })
      }

      for(j in 1:length(trees[[i]])){
        trees[[i]][[j]]@tree <- expanded_list[[j]]
      }

    }

  }else if(inherits(trees, 'AntibodyForests')){
    trees@tree <- trees %>% expand_intermediates()


  }else{
    stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
  }

  return(trees)
}
