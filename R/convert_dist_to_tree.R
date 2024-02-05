#' Converts string distance matrix into B cell lineage tree
#' Authors: Valentijn Tromp
#' @description This function takes a pairwise string distance matrix representing distances between B cell nodes and converts it into a lineage tree using various distance-based tree construction algorithms. 
#' @param clone string - specifies the sample and clonotype of which the distance matrix is derived. It should be in the format 'S1_clonotype1'.
#' @param dist.matrix matrix - contains string distances between all possible pairs of nodes. Row- and colnames correspond to the node names.
#' @param network.algorithm string - denotes distance-based tree construction method to be used to convert the distance matrix into a lineage tree. Options: 'network.tree', 'network.mst', and 'network.nj'. Defaults to 'network.tree'.
#' 'network.tree' : mst-like tree evolutionary network algorithm in which the germline node is positioned at the top of the tree, and nodes with the minimum distance to any existing node in the tree are linked iteratively.
#' 'network.mst'  : minimum spanning tree (MST) algorithm from 'ape::mst()' constructs networks/trees with the minimum sum of edge lengths/weights, which involves iteratively adding edges to the network in ascending order of edge weights, while ensuring that no cycles are formed.
#' 'phylo.nj'     : neighbor-joining (NJ) algorithm from 'ape::nj()' constructs phylogenetic trees by joining pairs of nodes with the minimum distance, creating a bifurcating tree consisting of internal nodes (representing unrecovered sequences) and terminal nodes (representing the recovered sequences). After tree construction, the tree is rerooted with respect to the germline node and internal nodes having a distance of 0 to a terminal node are replaced by that terminal node.
#' @param resolve.ties string or vector of strings - denotes the way ties are handled during the conversion of the distance matrix into lineage trees by the 'network.tree' algorithm (in the event where an unlinked node, that is to be linked to the tree next, shares identical distances with multiple previously linked nodes in the lineage tree). Options: 'min.expansion', 'max.expansion', 'min.germline.dist', 'max.germline.dist', 'min.germline.edges', 'max.germline.edges', and 'random'. If a vector is provided, ties will be resolved in a hierarchical manner. Defaults to 'c("max.expansion", "close.germline.dist", "close.germline.edges")'.
#' 'min.expansion'        : the node(s) having the smallest size is/are selected.
#' 'max.expansion'        : the node(s) having the biggest size is/are selected.
#' 'min.germline.dist'    : the node(s) having the smallets string distance to the germline node is/are selected.
#' 'max.germline.dist'    : the node(s) having the biggest string distance to the germline node is/are selected.
#' 'min.germline.edges'   : the node(s) having the lowest possible number of edges to the germline node is/are selected.
#' 'max.germline.edges    : the node(s) having the highest possible number of edges to the germline node is/are selected.
#' 'random'               : a random node is selected.
#' @param node.sizes list - contains node labels as indices and node sizes as corresponding values. These sizes are appended to the node attributes of igraph object and used by the 'min.expansion' and 'max.expansion' option to handle ties.
#' @return Returns a B cell lineage tree in the form of an igraph object. 
#' @export
#' @examples
#' \dontrun{
#' convert_dist_to_tree(dist.matrix = string_distance_matrix,
#'                      network.algorithm = "tree",
#'                      resolve.ties = c("max.expension", min.germline.dist", "min.germline.edges"),
#'                      node.sizes = node_sizes)
#'}


convert_dist_to_tree <- function(clone,
                                 dist.matrix, 
                                 network.algorithm = "network.tree", 
                                 resolve.ties = c("max.expansion", "close.germline.dist", "close.germline.edges"), 
                                 node.sizes){
  
  
  count_edges_between_two_nodes <- function(node1, node2, edge.matrix){
    
    # Counts edges to go from node1 to node2 from a matrix containing edges of a network/tree
    # Arguments:
    # - node1: the starting node (up in the network/tree)
    # - node2: the ending node (down in the network/tree)
    # - edge_matrix: 2-column matrix containing nodes of a directed network in which each row represent an edge, while the first column contains the upper node and the second column contains the lower node
    # Author: Valentijn Tromp
    
    # Save input edge matrix in 'edge_matrix'
    edge_matrix <- edge.matrix
    
    # Count number of edges by finding path from 'node2' to 'node1' (thereby moving upwards in the network)
    current_node <- node2
    
    # Initialize counter for the edges
    edge_count <- 0
    
    # Keep counting edges until the 'current_node' is equal to 'node1'
    while(current_node != node1){
      
      # Find node to which the 'current_node' is connected (the 'current_node' will be present in second column, while the node to which the 'current_node' is connected will be present in the first column)
      next_node <- edge_matrix[edge_matrix[, 2] == current_node, 1]
      
      # Update 'edge_count' and 'current_node'
      edge_count <- edge_count+1
      current_node <- next_node
    }
    
    # Return 'edge_count' 
    return(edge_count)
  }
  
  
  # 0. Make back-up of distance matrix and remove diagonal 0's 
  
  # Save distance matrix in 'dist_matrix'
  dist_matrix <- dist.matrix
  
  # Save node sizes in 'node_sizes'
  node_sizes <- node.sizes
  
  # Duplicate 'dist_matrix' for later on
  dist_matrix_backup <- dist_matrix
  
  
  # 1. Create lineage tree with 'network.tree' algorithm and store edges with edge lengths/weights in 'edges' dataframe (if 'network.algorithm' is set to 'network.tree')
  
  if(network.algorithm == "network.tree"){
    
    # Create matrix to store edges of the tree 
    edges <- matrix(nrow = 0, ncol = 3)
    colnames(edges) <- c("upper.node", "lower.node", "edge.length")
    
    # Replace diagonal 0's by NA
    diag(dist_matrix) <- NA
    
    # Get minimum distance to the germline node
    min_distance <- min(na.exclude(dist_matrix[,"germline"]))
    
    # Select nodes that have this distance to the germline node
    nodes_to_be_connected <- rownames(dist_matrix)[dist_matrix[, "germline"] == min_distance & !is.na(dist_matrix[, "germline"])]
    
    # Add edge(s) between germline and node(s) with 'min_distance' to the germline to the 'edges' matrix
    for(node in nodes_to_be_connected){
      edges <- rbind(edges, c("germline", node, dist_matrix_backup["germline", node]))
    }
    
    # Remove minimum distance to the germline node from the distance matrix
    dist_matrix[,"germline"][dist_matrix[,"germline"] == min_distance] <- NA
    dist_matrix["germline",][dist_matrix["germline",] == min_distance] <- NA
    
    # Iteratively add other nodes to the 'edges' matrix until all nodes are present in the 'edges' matrix
    while(all(rownames(dist_matrix) %in% edges) == FALSE){
      
      # Select nodes that are already connected (present in the 'edges' matrix)
      already_connected_nodes <- unique(na.exclude(c(edges[,1], edges[,2])))
      
      # Get minimum distance to one of these nodes
      min_distance <- min(na.exclude(dist_matrix[,already_connected_nodes]))
      
      # Select unlinked nodes that have this distance to any existing node in 'already_connected_nodes'
      nodes_to_be_connected <- unique(unlist(lapply(rownames(dist_matrix), function(x) if(min_distance %in% dist_matrix[x,already_connected_nodes] && !(x %in% already_connected_nodes)){return(x)})))
      
      # Iterate through nodes in 'nodes_to_be_connected'
      for(node in nodes_to_be_connected){
        
        # Select the node(s) in 'already_connected_nodes' that has 'min_distance' to the current node
        node_to_connect_to <- unlist(lapply(already_connected_nodes, function(x) if(dist_matrix[x, node] == min_distance){return(x)}))
        
        # If multiple nodes in 'node_to_connect_to' share 'min_distance' to the current node, the 'resolve.ties' options are hierarchically applied to narrow down the candidates for linking the currently unlinked node in the tree to, aiming to select a single node.
        if(length(node_to_connect_to) != 1){
          
          # Loop over provided 'resolve.ties' options
          for(option in resolve.ties){
            
            # If the 'resolve.ties' parameter is set to 'min.expansion', the node(s) having the smallest size is/are selected
            if(option == "min.expansion"){
              node_to_connect_to <- node_to_connect_to[node_sizes[node_to_connect_to] == min(unlist(node_sizes[node_to_connect_to]))]
            }
            
            # If the 'resolve.ties' parameter is set to 'max.expansion', the node(s) having the biggest size is/are selected
            if(option == "max.expansion"){
              node_to_connect_to <- node_to_connect_to[node_sizes[node_to_connect_to] == max(unlist(node_sizes[node_to_connect_to]))]
            }
            
            # If the 'resolve.ties' parameter is set to 'min.germline.dist', the node(s) having the smallest string distance to the germline node is/are selected
            if(option == "min.germline.dist"){
              node_to_connect_to <- node_to_connect_to[dist_matrix_backup[node_to_connect_to, "germline"] == min(dist_matrix_backup[node_to_connect_to, "germline"])]
            }
            
            # If the 'resolve.ties' parameter is set to 'max.germline.dist', the node(s) having the biggest string distance to the germline node is/are selected
            if(option == "max.germline.dist"){
              node_to_connect_to <- node_to_connect_to[dist_matrix_backup[node_to_connect_to, "germline"] == max(dist_matrix_backup[node_to_connect_to, "germline"])]
            }
            
            # If the 'resolve.ties' parameter is set to 'min.germline.edges', the node(s) having the lowest possible number of edges to the germline node is/are selected
            if(option == "min.germline.edges"){
              number_of_edges <- sapply(node_to_connect_to, function(node) count_edges_between_two_nodes(node1 = "germline", node2 = node, edge.matrix = edges))
              node_to_connect_to <- node_to_connect_to[number_of_edges == min(number_of_edges)]
            }
            
            # If the 'resolve.ties' parameter is set to 'max.germline.edges', the node(s) having the highest possible number of edges to the germline node is/are selected
            if(option == "max.germline.edges"){
              number_of_edges <- sapply(node_to_connect_to, function(node) count_edges_between_two_nodes(node1 = "germline", node2 = node, edge.matrix = edges))
              node_to_connect_to <- node_to_connect_to[number_of_edges == max(number_of_edges)]
            }
            
            # If the 'resolve.ties' parameter is set to 'random', a random node will be selected
            if(option == "random"){
              node_to_connect_to <- sample(node_to_connect_to, 1)
            }
          }
        }
        
        # If there is only one node left in 'node_to_connect_to'
        if(length(node_to_connect_to) == 1){
          
          # Add new edge between this linked node in the tree and the currently unlinked node to the 'edges' matrix
          edges <- rbind(edges, c(node_to_connect_to, node, dist_matrix_backup[node_to_connect_to, node]))
        }
        
        # If there are still multiple nodes in 'node_to_connect_to' that has 'min_distance' to the current node...
        if(length(node_to_connect_to) != 1){
          
          # Print warning message
          warning(paste0(c("Not all ties could be resolved in ", clone, " !!!")))
          
          # The currently unlinked node is connected to both linked nodes in the tree, thereby creating a loop in the network
          for(i in node_to_connect_to){
            edges <- rbind(edges, c(i, node, dist_matrix_backup[i, node]))
          }
        }
        
        # Remove this minimum distance from the distance matrix
        dist_matrix[node_to_connect_to, node] <- NA
        dist_matrix[node, node_to_connect_to] <- NA
      }
    }
    
    # Convert 'edges' matrix into dataframe
    edges <- as.data.frame(edges)
    
  }
  
  
  # 2. Create lineage tree with 'mst' algorithm and store edges with edge lengths/weights in 'edges' dataframe (if 'network.algorithm' is set to 'mst')
  
  if(network.algorithm == "network.mst"){
    
    # Create minimum spanning stree with 'ape::mst()' function
    mst_network <- ape::mst(dist_matrix)
    
    # Convert 'mst' object into 'igraph' object
    mst_igraph_object <- igraph::graph_from_adjacency_matrix(mst_network, mode = "undirected")
    
    # Create matrix 'edges' containing the edges of the minimum spanning tree
    edges <- igraph::get.edgelist(mst_igraph_object)
    
    # Convert 'edges' matrix into dataframe and rename columns
    edges <- as.data.frame(edges)
    colnames(edges) <- c("upper.node", "lower.node")
    
    # Add a third column 'edge.lengths' containing the length/weight of the edges
    edges$edge.length <- sapply(1:nrow(edges), function(x) dist_matrix[edges[x, "upper.node"], edges[x, "lower.node"]])
    
    # Start reordering the nodes in the 'edges' dataframe from the 'germline'
    current_upper_nodes <- "germline"
    processed_nodes <- current_upper_nodes
    
    # Keep reordering the nodes in the 'edges' dataframe until all nodes are processed and in the 'processed_nodes' vector
    while(all(rownames(dist_matrix) %in% processed_nodes) == FALSE){
      
      # Initialize vector to store nodes that are added in this step
      processed_lower_nodes <- c()
      
      # Iterate through nodes in 'current_upper_nodes'
      for(upper_node in current_upper_nodes){
        
        # Select the rows/edges from 'edges' dataframe that contain the current 'upper_node'
        selected_edges <- edges[edges$upper.node == upper_node | edges$lower.node == upper_node, ]
        
        # Retrieve all the nodes that are present in this selection of edges
        selected_nodes <- unique(c(selected_edges$upper.node, selected_edges$lower.node))
        
        # Remove the nodes that are already processed, the remaining nodes will be the lower nodes of the current 'upper_node'
        current_lower_nodes <- selected_nodes[!selected_nodes %in% processed_nodes]
        
        # Iterate through nodes in 'current_lower_nodes'
        for(lower_node in current_lower_nodes){
          
          # Swap the nodes in the rows in the 'edges' dataframe. in which the current 'upper_node' is present in the 'lower.node' column and the current 'lower_node' is present in the 'upper.node' column
          edges[edges$lower.node == upper_node & edges$upper.node == lower_node, ] <- edges[edges$lower.node == upper_node & edges$upper.node == lower_node, c("lower.node", "upper.node", "edge.length")]
        }
        
        # After iterating trough nodes in 'current_lower_nodes' and swapping nodes in 'edges' dataframe, append 'current_lower_nodes' to 'processed_lower_nodes'
        processed_lower_nodes <- c(processed_lower_nodes, current_lower_nodes)
      }
      
      # All the nodes in 'processed_lower_nodes' will be the 'current_upper_nodes' in the next iteration
      current_upper_nodes <- processed_lower_nodes
      
      # Update 'processed_nodes' vector by appending nodes in 'processed_lower_nodes'
      processed_nodes <- c(processed_nodes, processed_lower_nodes)
    }
  }
  
  
  # 3. Create lineage tree with 'nj' algorithm and store edges with edge lengths/weights in 'edges' dataframe (if 'network.algorithm' is set to 'nj')
  
  if(network.algorithm == "phylo.nj"){
    
    # A neighbor joining tree can only be build if there are 3 or more sequences
    if(nrow(dist_matrix) > 2){
      
      # Create neighbor joining tree with 'ape::nj()' function
      nj_phylo_tree <- ape::nj(dist_matrix)
      
      # Reroot tree with respect to germline node
      nj_phylo_tree_rerooted <- ape::root(nj_phylo_tree, "germline")
      
      # Convert 'phylo' object into 'igraph' object
      nj_igraph_object <- ape::as.igraph.phylo(nj_phylo_tree_rerooted)
      
      # Create matrix containing the edges of neighbor joining tree between internal nodes (first column) and terminal nodes (second column)
      edges <- igraph::get.edgelist(nj_igraph_object)
      
      # Remove 'Node' from labels of internal nodes
      edges <- gsub(pattern = "Node", replacement = "", x = edges)
      
      # Create dataframe containing the edges plus a third column containing the length of the edges
      edges <- data.frame(do.call(rbind, lapply(1:nrow(edges), function(x) c(edges[x, 1], edges[x, 2], as.numeric(nj_phylo_tree_rerooted$edge.length[x])))))
      
      # Rename column names of 'edges' dataframe
      colnames(edges) <- c("upper.node", "lower.node", "edge.length")
      
      # Make a subset of edges with a length of 0
      zero_length_edges <- edges[edges$edge.length == 0,]
      
      # Remove edges with a length of 0 from 'edges' dataframe
      edges <- edges[edges$edge.length != 0,]
      
      # If the 'zero_length_edges' is not empty...
      if(nrow(zero_length_edges) != 0){
        
        # Loop over edges in 'zero_length_edges' dataframe
        for(i in 1:nrow(zero_length_edges)){
          
          # Retrieve internal and terminal node that are connected by this 'zero length edge'
          internal_node <- zero_length_edges[i, "upper.node"]
          terminal_node <- zero_length_edges[i, "lower.node"]
          
          # Replace 'internal_node' with 'terminal_node' in the columns 'upper.node' and 'lower.node' in the 'edges' dataframe
          edges[edges$upper.node == internal_node, "upper.node"] <- terminal_node
          edges[edges$lower.node == internal_node, "lower.node"] <- terminal_node
        }
      }
      
      # Swap values in rows where the second column contains the 'germline' node (in order to place the germline node on top of the tree)
      edges[edges$lower.node == "germline", ] <- edges[edges$lower.node == "germline", c("lower.node", "upper.node", "edge.length")]
    }
    
    # If there are less than 3 sequences, the tree consists of a germline node and a single descendant 
    if(nrow(dist_matrix) < 3){
      edges <- data.frame(upper.node = "germline", lower.node = "node1", edge.length = dist_matrix["germline", "node1"])
    }
  }
  
  
  # 4. Convert 'edges' dataframe into igraph object  
  
  # Create 2-column matrix in which each row represent an edge in the graph
  edge_matrix <- as.matrix(edges[, c("upper.node", "lower.node")])
  
  # Create igraph object using edges stored in 'edge_matrix'
  tree <- igraph::graph_from_edgelist(edge_matrix)
  
  # Append node size from 'node_sizes' ligt to vertex/node attributes of igraph object
  igraph::V(tree)$size <- node_sizes[igraph::V(tree)$name]
  
  # Append edge length/weight from 'edges' dataframe to edge attributes of igraph object
  igraph::E(tree)$length <- edges$edge.length
  
  # Return tree
  return(tree)
}