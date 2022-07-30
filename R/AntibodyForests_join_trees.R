#' Joins a list of trees/networks as AntibodyForests objects into a single AntibodyForests object

#'@description Joins a list of trees/networks as AntibodyForests into a single AntibodyForests object. The resulting network will include all joined networks as separate components. Useful for faster downstream analyses (e.g., node metrics via AntibodyForests_metrics on this object instead on each separate object, plotting multiple trees in the same plot, etc.,)
#' @param tree.list (nested) list of AntibodyForests objects, as obtained from the AntibodyForests function.
#' @param join.per string - 'sample' joins the objects per sample if the input is a nested list of AntibodyForests objects, resulting in a single joined graph per sample, 'global' joins all graphs in the nested list into a single object.
#' @param join.method string - networks, especially minimum spanning trees, can be joined into a single connected graph if join.method = 'single.germline' (will pick a single germline from all available germlines and will recalculate the string distance) or 'multiple.germlines.joined' (will add inter-germline edges).
#' Will create a single object with disconnected subgraph if join.method = 'multiple.germlines'.

#' @return single AntibodyForests object consisting of the joined graphs/trees. Resulting graph can be either a single connected component (if join.method = 'single.germline' or 'multiple.germlines.joined') or multiple disconnected subgraphs (join.method = 'multiple.germlines').
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_join_trees(tree.list, join.per = 'sample', join.method = 'multiple.germlines')
#'}

AntibodyForests_join_trees <- function(tree.list, join.per, join.method){

  if(missing(tree.list)) stop('Please input a (nested) list of AntibodyForests objects which you wish to join.')
  if(missing(join.per)) join.per <- 'sample'
  if(missing(join.method)) join.method <- 'multiple.germlines'

  new <- NULL


  methods::setClass('AntibodyForests',
    slots = c(
      tree = 'ANY', #in
      sample_id = 'ANY', #in
      clonotype_id = 'ANY', #in
      plot_ready = 'ANY', #f call
      heterogeneous = 'ANY', #f call
      reactivity = 'ANY', #f call
      dynamic = 'ANY', #f call
      metrics = 'ANY', #f call
      sequences = 'ANY', #in
      germline_sequence = 'ANY', #to add
      barcodes = 'ANY', #in
      node_features = 'ANY', #in
      edge_list = 'ANY', #no/f call
      gex_list = 'ANY', #f call
      paths = 'ANY', #f call
      node_transitions = 'ANY', #f call
      network_algorithm = 'ANY', #in
      adjacency_matrix = 'ANY', #no/f call
      phylo = 'ANY',
      feature_names = 'ANY',
      permuted_transitions = 'ANY' #f call
    )
  )


  join_trees <- function(tree.list){
    directed <- igraph::is_directed(tree.list[[1]]@tree)

    adjacency_matrix_list <- lapply(tree.list, function(x) as.matrix(igraph::as_adjacency_matrix(x@tree, attr = 'weight')))
    vertex_dfs_list <- lapply(tree.list, function(x) igraph::as_data_frame(x@tree, what = 'vertices'))

    if(join.method=='multiple.germlines' | join.method=='multiple.germlines.joined'){

      output_df <- do.call('rbind', vertex_dfs_list)
      nodes_per_network <- unname(lapply(adjacency_matrix_list, function(x) nrow(x)))
      output_matrix_nodes <- sum(unlist(nodes_per_network))
      output_matrix <- matrix(0, output_matrix_nodes, output_matrix_nodes)
      #output_df$sequence_id <- 1:nrow(output_df)

      germline_nodes <- which(output_df$node_type=='germline')
      m1 <- adjacency_matrix_list[[1]]
      prev_added_node <- 1
      current_node <- nrow(m1)
      output_matrix[prev_added_node:current_node, prev_added_node:current_node] <- m1[1:nrow(m1), 1:nrow(m1)]
      prev_added_node <- current_node + 1
      m1 <- NULL

      for(m in adjacency_matrix_list[-1]){
        current_node <- prev_added_node + nrow(m) - 1
        output_matrix[prev_added_node:current_node, prev_added_node:current_node] <- m
        prev_added_node <- current_node + 1
        m<-NULL
      }

      if(join.method=='multiple.germlines.joined'){
        germline_nodes <- which(output_df$node_type=='germline')
        germlines <- output_df$network_sequences[germline_nodes]

        distance_matrix <- stringdist::stringdistmatrix(germlines, germlines, method='lv')

        for(i in 1:(length(germline_nodes)-1)){
          edge_weight <- distance_matrix[i,i+1]


          output_matrix[germline_nodes[i], germline_nodes[i+1]] <- edge_weight
          if(!directed){
            output_matrix[germline_nodes[i+1], germline_nodes[i]] <- edge_weight
          }
        }
      }

    }else if(join.method=='single.germline'){
      output_df <- do.call('rbind', vertex_dfs_list)

      nodes_per_network <- lapply(adjacency_matrix_list, function(x) nrow(x))
      germlines_per_network <- lapply(vertex_dfs_list, function(x) any(x$node_type=='germline'))

      output_matrix_nodes <- sum(unlist(nodes_per_network)) - length(which(germlines_per_network==T)) + 1
      output_matrix <- matrix(0, output_matrix_nodes, output_matrix_nodes)

      germline_nodes_new <- unlist(which(output_df$node_type=='germline'))
      germline_nodes_original <- as.integer(output_df$label[germline_nodes_new])

      germline_row <- output_df[germline_nodes_new[1], ]
      rownames(germline_row) <- NULL
      output_df <- output_df[which(output_df$node_type !='germline'),]
      rownames(output_df) <- NULL
      output_df <- rbind(output_df, germline_row)

      new_germline_node <- nrow(output_matrix)
      new_germline_sequence <- output_df$network_sequences[new_germline_node]

      nodes_connected_to_germline <- mapply(function(x,y) which(x[y,]!=0), adjacency_matrix_list, germline_nodes_original)
      sequences_connected_to_germline <- mapply(function(x,y) x$network_sequences[y], vertex_dfs_list, nodes_connected_to_germline)
      new_edge_weights <- lapply(sequences_connected_to_germline, function(x) stringdist::stringdistmatrix(unlist(x), new_germline_sequence))

      nodes_connected_to_germline
      m1 <- adjacency_matrix_list[[1]]
      prev_added_node <- 1
      current_node <- nrow(m1)
      current_node <- current_node - 1
      output_matrix[prev_added_node:current_node, prev_added_node:current_node] <- m1[1:(nrow(m1)-1), 1:(nrow(m1)-1)]
      prev_added_node <- current_node + 1
      m1 <- NULL
      matrices <- adjacency_matrix_list[-1]
      output_matrix[new_germline_node,nodes_connected_to_germline[[1]]] <- new_edge_weights[[1]]

      if(!directed){
        output_matrix[nodes_connected_to_germline[[1]], new_germline_node] <- new_edge_weights[[1]]
      }
      #matrices[[i]][1:(nrow(matrices[[i]])-1), 1:(nrow(matrices[[i]])-1)]
      for(i in 1:length(matrices)){
        current_node <- prev_added_node + nrow(matrices[[i]]) - 1
        current_node <- current_node - 1
        output_matrix[prev_added_node:current_node, prev_added_node:current_node] <-  matrices[[i]][1:(nrow(matrices[[i]])-1), 1:(nrow(matrices[[i]])-1)]

        output_matrix[new_germline_node, unlist(nodes_connected_to_germline[i+1]) + prev_added_node - 1] <- new_edge_weights[[i+1]]
        if(!directed){
          output_matrix[unlist(nodes_connected_to_germline[i+1]) + prev_added_node - 1, new_germline_node] <- new_edge_weights[[i+1]]
        }
        prev_added_node <- current_node + 1
      }
    }

    edgelist  <-  igraph::graph_from_adjacency_matrix(output_matrix, weighted = TRUE) %>%
                  igraph::as_data_frame()


    output_df$label <- 1:nrow(output_df)
    output_df$name <- 1:nrow(output_df)

    g <- igraph::graph_from_data_frame(edgelist, directed = directed, vertices = output_df)
    igraph::V(g)$label <- 1:length(igraph::V(g))

    return(g)
  }


  instantiate_tree_object <- function(g){

    node_df <- igraph::as_data_frame(g, what = 'vertices')
    node_df_subset <- node_df[node_df$node_type == 'sequence',]
    #edge_df <- igraph::as_data_frame(g, what = 'edges')

    labels <- node_df$labels

    network_sequences <- node_df_subset$network_sequences
    names(network_sequences) <- labels

    cell_barcodes <- node_df_subset$cell_barcodes
    names(cell_barcodes) <- labels

    sample_id <- unlist(unique(node_df$sample_id))
    sample_id <- sample_id[order(nchar(sample_id), sample_id)]
    clonotype_id <- unlist(unique(node_df$clonotype_id))
    clonotype_id <- clonotype_id[order(nchar(clonotype_id), clonotype_id)]

    ####optional
    node_df$network_sequences <- NULL
    node_df$cell_barcodes <- NULL
    node_df$sample_id <- NULL
    node_df$clonotype_id <- NULL
    #node_df$node_type <- NULL

    list_column_names <- names(sapply(node_df, class)[which(sapply(node_df, class)=='list')])
    list_column_names <- list_column_names[stringr::str_detect(list_column_names, '_counts')]
    list_column_names <- stringr::str_remove(list_column_names, pattern='_counts')

    for(col in list_column_names){
      node_df <- node_df[,!(names(node_df) %in% paste0(col, '_counts'))]
    }


    final_object <- new('AntibodyForests',
                        tree = g,
                        sample_id = sample_id,
                        clonotype_id = clonotype_id,
                        sequences = network_sequences,
                        barcodes = cell_barcodes,
                        node_features = node_df,
                        feature_names = feature.names,
                        network_algorithm = network.algorithm
                        #edge_list = edge_df,
                        #adjacency_matrix = igraph::as_adjacency_matrix(g),
                       )

    return(final_object)
  }


  if(inherits(tree.list[[1]], 'AntibodyForests')){
    feature.names <- tree.list[[1]]@feature_names
    network.algorithm <- tree.list[[1]]@network_algorithm

    final_object <- tree.list %>%
                    join_trees() %>%
                    instantiate_tree_object()


  }else if(join.per == 'sample'){

    final_object <- vector(mode = "list", length = length(tree.list))

    for(i in 1:length(tree.list)){
      feature.names <- tree.list[[i]][[1]]@feature_names
      network.algorithm <- tree.list[[i]][[1]]@network_algorithm

      final_object[[i]] <- tree.list[[i]] %>%
                           join_trees() %>%
                           instantiate_tree_object()
    }

    if(length(final_object) == 1){
      final_object <- final_object[[1]]
    }


  }else if(join.per == 'global'){
    feature.names <- tree.list[[1]][[1]]@feature_names
    network.algorithm <- tree.list[[1]][[1]]@network_algorithm

    final_object <- purrr::flatten(tree.list) %>%
                         join_trees() %>%
                         instantiate_tree_object()


  }else{
    stop('Unrecognized joining method!')
  }

  return(final_object)

}
