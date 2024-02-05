#' Infers B cell evolutionary networks for all clonotypes in VDJ dataframe as obtained from the minimal_VDJ() function.
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description The resulting B cell lineage trees/networks provide insights into the evolutionary relationships between B cell sequences from each clonotype.
#' @param VDJ dataframe - VDJ object as obtained from the minimal_VDJ() function in Platypus.
#' @param sequence.columns string or vector of strings - denotes the sequence column(s) in the VDJ dataframe that contain the sequences that will be used to infer B cell lineage trees. Nodes in the trees will represent unique combinations of the selected sequences.
#' @param germline.columns string or vector of strings - denotes the germline column(s) in the VDJ dataframe that contain the sequences that will be used as starting points of the trees. The columns should be in the same order as in 'sequence.columns'.
#' @param concatenate.sequences bool - if TRUE, sequences from multiple sequence columns are concatenated into one sequence for single distance matrix calculations / multiple sequence alignments, else, a distance matrix is calculated / multiple sequence alignment is performed for each sequence column separately. Defaults to FALSE. 
#' @param node.features string or vector of strings - denotes the column name(s) in the VDJ dataframe from which the node features should be extracted (which can, for example, be used for plotting of lineage trees later on).
#' @param network.algorithm string - denotes the network algorithm that will be used to convert the distance matrix or multiple sequence alignment into a lineage tree/network. Distance-based options: 'network.tree', 'network.mst', and 'phylo.nj'. Multiple sequence alignment-based options: 'phylo.ml', and 'phylo.mp'. Defaults to 'network.tree' (mst-like algorithm).
#' 'network.tree' : mst-like tree evolutionary network algorithm in which the germline node is positioned at the top of the tree, and nodes with the minimum distance to any existing node in the tree are linked iteratively.
#' 'network.mst'  : minimum spanning tree (MST) algorithm from 'ape::mst()' constructs networks/trees with the minimum sum of edge lengths/weights, which involves iteratively adding edges to the network in ascending order of edge weights, while ensuring that no cycles are formed.
#' 'phylo.nj'     : neighbor-joining (NJ) algorithm from 'ape::nj()' constructs phylogenetic trees by joining pairs of nodes with the minimum distance, creating a bifurcating tree consisting of internal nodes (representing unrecovered sequences) and terminal nodes (representing the recovered sequences). After tree construction, the tree is rerooted with respect to the germline node and internal nodes having a distance of 0 to a terminal node are replaced by that terminal node.
#' 'phylo.ml'     : --- TO BE ADDED SOON! ---
#' 'phylo.mp'     : --- TO BE ADDED SOON! ---
#' @param distance.metric string - denotes the metric that will be calculated with the 'stringdist()' function to measure (string) distance between sequences. Options: 'lv', 'dl', 'osa', 'hamming', 'lcs', 'qgram', 'cosine', 'jaccard', 'jw', and 'soundex'.  Defaults to 'lv' (Levenshtein distance / edit distance).
#' 'lv'       : Levensthein distance (also known as edit distance) equals to the minimum number of single-element edits (insertions, deletions, or substitutions) required to transformer one string into another.
#' 'dl'       : Damerau-Levenshtein distance is similar to the Levenshtein distance, but also allows transpositions of adjacent elements as a single-edit operation.
#' 'osa'      : Optimal String Alignment distance is similar to the Damerau-Levensthein distance, but does not allow to apply multiple transformations on a same substring.
#' 'hamming'  : Hamming distance equals to the number of positions at which the corresponding elements differ between two strings (applicable only to strings of equal length).
#' 'lcs'      : Longest Common Subsequence distance is similar to the Levenshtein distance, but only allowing insertions and deletions as single-edit operations.
#' 'qgram'    : Q-gram distance equal to the number of distinct q-grams that appear in either string but not both, whereby q-grams are all possible substrings of length q in both strings (q defaults to 1).
#' 'cosine'   : cosine distance equals to 1 - cosine similarity (the strings are converted into vectors containing the frequency of all single elements (A and B), whereby the cosine similarity (Sc) equals to the dot product of these vectors divided by the product of the magnitude of these vectors, which can be written in a formula as Sc(A, B) = A . B / (||A|| x ||B||)).
#' 'jaccard'  : Jaccard distance equals to 1 - Jaccard index (the strings are converted into sets of single elements (A and B), whereby the Jaccard index (J) equals to the size of the intersection of the two sets divided by the size of the union of the sets, which can be written in a formula as J(A, B) = |A ∩ B| / |A ∪ B|).
#' 'jw'       : Jaro-Winkler distance equals to 1 - Jaro-Winkler similarity (Jaro-Winkler similary is calculated with the following formulas: Sw = Sj + P * L * (1-Sj) in which Sw is the Jaro-Winkler similary, Sj is the Jaro similarity, P is the scaling factor (defaults to 0), and L is the length of the matching prefix; and Sj = 1/3 * (m/|s1| + m/|s2| + (m-t)/m) in which Sj is the Jaro similarity, m is the number of matching elements, |s1| and |s2|are the lengths of the strings, and t is the number of transpositions).
#' 'soundex'  : Soundex distance equals to 1 if the 4-character Soundex code of the strings do not match (Soundex is a phonetic algorithm that converts strings into a 4-character code based on their (English) pronunciation).
#' @param substitution.model string - denotes the codon substitution model that will be used during the likelihood calculations (if 'network.algorithm' is set to 'phylo.ml'. Options: 'GY94' and 'HLP17'. Defaults to NULL. 
#' 'GY94'   : --- TO BE ADDED SOON! ---
#' 'HLP17'  : --- TO BE ADDED SOON! ---
#' @param resolve.ties string or vector of strings - denotes the way ties are handled during the conversion of the distance matrix into lineage trees by the 'network.tree' algorithm (in the event where an unlinked node, that is to be linked to the tree next, shares identical distances with multiple previously linked nodes in the lineage tree). Options: 'min.expansion', 'max.expansion', 'min.germline.dist', 'max.germline.dist', 'min.germline.edges', 'max.germline.edges', and 'random'. If a vector is provided, ties will be resolved in a hierarchical manner. Defaults to 'c("max.expansion", "close.germline.dist", "close.germline.edges")'.
#' 'min.expansion'        : the node(s) having the smallest size is/are selected.
#' 'max.expansion'        : the node(s) having the biggest size is/are selected.
#' 'min.germline.dist'    : the node(s) having the smallets string distance to the germline node is/are selected.
#' 'max.germline.dist'    : the node(s) having the biggest string distance to the germline node is/are selected.
#' 'min.germline.edges'   : the node(s) having the lowest possible number of edges to the germline node is/are selected.
#' 'max.germline.edges    : the node(s) having the highest possible number of edges to the germline node is/are selected.
#' 'random'               : a random node is selected.
#' @param parallel bool - if TRUE, the per-clone network inference is executed in parallel (parallelized across samples). Defaults to FALSE.
#' @param num.cores integer - number of cores to be used when parallel = TRUE. Defaults to all available cores - 1 or the number of samples in the VDJ dataframe (depending which number is smaller).
#' @return Returns nested list of AntibodyForests objects for each sample and each clonotype, whereby an AntibodyForests object consists of a list of all nodes with their sequences and features called 'clones', a list of distance matrices, one for each selected sequence column, called 'distance.matrices', and the lineage tree as an igraph object called 'lineage.tree'. For example, output[[1]][[2]] denotes the AntibodyForests object of the first sample, second clonotype. 
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests(VDJ, 
#'                 sequence.columns = c("VDJ_sequence_aa_trimmed","VJ_sequence_aa_trimmed"),
#'                 germline.columns = c("VDJ_germline_aa_trimmed","VJ_germline_aa_trimmed"), 
#'                 concatenate.sequences = TRUE,
#'                 node.features = c("VDJ_vgene", "VDJ_dgene", "VDJ_jgene", "VDJ_cgene"),
#'                 distance.metric = "lv",
#'                 network.algorithm = "tree",
#'                 resolve.ties = c("max.expansion", "close.germline.dist", "close.germline.edges"),
#'                 parallel = TRUE)
#'}


AntibodyForests <- function(VDJ,
                            sequence.columns,
                            germline.columns,
                            concatenate.sequences = FALSE,
                            node.features,
                            distance.metric,
                            substitution.model,
                            network.algorithm,
                            resolve.ties,
                            parallel = FALSE,
                            num.cores){
  
  # If 'sequence.columns' is missing or if the columns are not all present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(missing(sequence.columns) | !all(sequence.columns %in% colnames(VDJ))){stop("Error: Please provide valid 'sequence.columns' from the input VDJ dataframe.")}
  # If 'germline.columns' is missing or if the columns are not all present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(missing(germline.columns) | !all(germline.columns %in% colnames(VDJ))){stop("Error: Please provide valid 'germline.columns' from the input VDJ dataframe.")}
  # If the columns in 'germline.columns' do not correspond to the columns in 'sequence.columns', a warning is returned.
  if(!all(sequence.columns == gsub("germline", "sequence", germline.columns))){warning("WARNING: Please make sure that the columns in 'germline.columns' correspond to the columns in 'sequence.columns'.")}
  # If 'node.features' contains columns that are not present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(!all(node.features %in% colnames(VDJ))){stop("Error: Not all columns in 'node.features' could be found in the input VDJ dataframe.")}
  # If 'network.algorithm', 'distance.metric', 'substitution.model', and 'resolve.ties' are all missing, the edit distance / Levenshtein distance will be used to calculate pairwise string distance matrices and the 'network.tree' algorithm will be used to convert the distance matrices into lineage trees using the 'max.expansion', 'min.germline.dist', and 'min.germline.edges' options to resolve ties
  if(missing(network.algorithm) && missing(distance.metric) && missing(substitution.model) && missing(resolve.ties)){network.algorithm <- "network.tree"; distance.metric <- "lv"; resolve.ties <- c("max.expansion", "min.germline.dist", "min.germline.edges")}
  # If 'network.algorithm' is set to 'network.tree', but 'distance.metric' is missing
  if(network.algorithm == "network.tree" && missing(distance.metric)){distance.metric <- "lv"}
  # If 'network.algorithm' is set to 'network.tree', but 'resolve.ties' is missing, the 'max.expansion', 'min.germline.dist', and 'min.germline.edges' options are used to handle ties
  if(network.algorithm == "network.tree" && missing(resolve.ties)){resolve.ties <- c("max.expansion", "min.germline.dist", "min.germline.edges")}
  # If 'parallel' is set to TRUE but 'num.cores' is not specified, the number of cores is set to all available cores - 1
  if(parallel == TRUE && missing(num.cores)){num.cores <- parallel::detectCores() -1}
  
  
  convert_dist_to_tree <- function(clone,
                                   dist.matrix, 
                                   network.algorithm, 
                                   resolve.ties, 
                                   node.sizes){
    
    # Converts a pairwise distance matrix into a lineage tree 
    # Arguments:
    # - clone: string specifying the sample and clonotype in the format "S1_clonotype1"
    # - dist.matrix: matrix containing distances between all possible pairs of sequences/nodes
    # - network.algorithm: string specifying distance-based tree construction method
    # - resolve.ties: string or vector of strings denoting the way ties are handled in the 'network.tree' algorithm
    # - node.sizes: list containing the sequence/node sizes/frequencies
    # Author: Valentijn Tromp
    
    
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
  
  
  infer_single_network <- function(VDJ,
                                   clone,
                                   sequence.columns,
                                   germline.columns,
                                   concatenate.sequences,
                                   node.features,
                                   network.algorithm,
                                   distance.metric,
                                   substitution.model,
                                   resolve.ties){
    
    # Infers network/tree for a single clone within one sample
    # Arguments:
    # - VDJ: VDJ dataframe as obtained from the minimal_VDJ() function in Platypus
    # - clone: string denoting the sample ID and the clonotype ID for which the network will be inferred (in the format of "S1_clonotype4" or "S2_clonotype3")
    # - sequence.columns: string or vector of strings denoting the sequence columns in the 'VDJ' dataframe that contain the sequences that will be used to infer B cell lineage trees
    # - germline.columns: string or vector of strings denoting the germline columns in the 'VDJ' dataframe that contain the sequences that will be used as starting points of the lineage trees
    # - concatenate.sequences: bool indicating whether sequences from different sequence columns should be concatenated before pairwise distance calculations or multiple string alignment
    # - node.features: string or vector of strings denoting the additional column name(s) in the VDJ dataframe to be selected
    # - network.algorithm: string denoting the network algorithm that will be used to convert the distance matrices or multiple sequence alignments into lineage trees
    # - distance.metric: string denoting the metric that will be calculated to measure (string) distances between sequences (if no substitution model is specified)
    # - substitution.model: denotes codon substitution model that will be used to calculate the alignment score when comparing nucleotide sequences
    # - resolve.ties: string denoting the way ties are handled when the 'tree' network algorithm is being used 
    
    # 1. Create a nested list 'nodes_list', in which each sublist represents a clone/node in the graph, and each node represent a unique combination of sequences of the columns selected in 'sequence.columns'
    # Each node/sublist contains the following items:
    # - a (separate) string for all columns selected in 'sequence.columns' containing the sequence for that node
    # - barcodes: a list of all barcodes of the cells that have the unique combination of sequences
    # - frequency: an integer that corresponds to the number of cells (equal to the length of the list 'barcodes')
    # - a (separate) list for all columns selected in 'node.features'
    
    # Retrieve sample ID and clonotype ID from 'clone'
    sample <- strsplit(clone, split="_")[[1]][1]
    clonotype <- strsplit(clone, split="_")[[1]][2]
    
    # Store 'sequence.columns' and 'germline.columns' in vectors 'sequence_columns' and 'germline_columns', respectively
    sequence_columns <- sequence.columns
    germline_columns <- germline.columns
    
    # Make a subset of the VDJ dataframe with all cells from the specified sample and clone and only the columns 'barcode', the specified 'sequence.columns' and 'germline.columns', and additional columns specified in 'node.features'
    vdj_subset <- base::subset(VDJ, sample_id == sample & clonotype_id == clonotype)[,c("barcode", sequence_columns, germline_columns, node.features)]
    
    # Create vectors specifying
    nt.columns <- colnames(vdj_subset)[stringr::str_detect(colnames(vdj_subset), "_nt")]
    aa.columns <- colnames(vdj_subset)[stringr::str_detect(colnames(vdj_subset), "_aa")]
    
    # Trim of last one or two nucleotides to make the length of the nucleotide sequences in 'vdj_subset' multiple of 3
    vdj_subset[nt.columns] <- lapply(vdj_subset[nt.columns], function(x){
      ifelse(nchar(x) %% 3 != 0, substring(x, 1, nchar(x) - (nchar(x) %% 3)))
    })
    
    # If 'concatenate.sequences' is set to TRUE...
    if(concatenate.sequences){
      
      # Concatenate sequences in 'sequence_columns' and 'germline_columns' and store concatenated sequences in 'concatenated_sequence' and 'concatenated_germline', respectively
      vdj_subset$concatenated_sequence <- sapply(1:nrow(vdj_subset), function(x) do.call(paste0, as.list(vdj_subset[x, sequence_columns])))
      vdj_subset$concatenated_germline <- sapply(1:nrow(vdj_subset), function(x) do.call(paste0, as.list(vdj_subset[x, germline_columns])))
      
      # Update column names in 'sequence_columns' and 'germline_columns' vectors
      sequence_columns <- "concatenated_sequence"
      germline_columns <- "concatenated_germline"
    }
    
    # Create a dataframe containing unique combinations of the sequences in 'sequence_columns'
    intraclonal_sequences_df <- as.data.frame(unique(vdj_subset[,sequence_columns]))
    colnames(intraclonal_sequences_df) <- sequence_columns
    
    # Create a list of sublists where each sublist represents a unique combination of sequences (and each unique combination of sequences will represent a node in the graph later on)
    nodes_list <- lapply(1:nrow(intraclonal_sequences_df), function(x){
      # Extract sequences for the current row and save in 'node_sublist'
      node_sublist <- lapply(sequence_columns, function(y) intraclonal_sequences_df[x,y])
      # Rename indices of 'node_sublist' to 'sequence_columns'
      names(node_sublist) <- sequence_columns
      # Find matching rows in 'vdj_subset' based on the current combination of sequences
      matching_rows <- sapply(1:nrow(vdj_subset), function(y) all(intraclonal_sequences_df[x,] == vdj_subset[y,sequence_columns])) 
      # Add barcodes and size/frequency information to 'node_sublist'
      node_sublist["barcodes"] <- list(vdj_subset$barcode[matching_rows])
      node_sublist["size"] <- sum(matching_rows)
      # Add 'node.features' to 'node_sublist' as well
      for(i in node.features){node_sublist[i] <- list(vdj_subset[,i][matching_rows])}
      # Return list of sublists 
      return(node_sublist)
    })
    
    # Order the list based on the size/frequency of the (sub)clones/nodes in descending order
    nodes_list <- nodes_list[order(sapply(nodes_list, function(x) x$size), decreasing = TRUE)]
    
    # Rename the sublists to 'node1', 'node2', etc.
    names(nodes_list) <- paste0("node", seq_along(nodes_list))
    
    # Create list of size/frequencies 
    node_sizes <- lapply(names(nodes_list), function(node) nodes_list[[node]][["size"]])
    names(node_sizes) <- names(nodes_list)
    node_sizes["germline"] <- 0
    
    # Add a 'germline' node to the list, containing the the most abundant germline sequences from the specified 'germline_columns' (most of the time, each germline column contains only one unique sequence)
    nodes_list$germline <- lapply(germline_columns, function(x) names(table(vdj_subset[,x]))[which.max(table(vdj_subset[,x]))])
    names(nodes_list$germline) <- sequence_columns
    
    
    # 2.1 Calculate a pairwise string distance matrix for each column in 'sequence_columns' using the specified 'distance.metric', if 'network.algorithm' is set to a distance-based tree/network construction method 
    
    # If 'network.algorithm' is set to 'network.tree', 'network.mst', or 'phylo.nj' and if the provided 'distance.metric' is one of the methods for distance calculation of the 'stringdistmatrix()' function, this function is used to calculate the distance matrix with the specified string distance metric
    if(network.algorithm %in% c("network.tree", "network.mst", "phylo.nj") && distance.metric %in% c("osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex")){
      
      # Calculate a distance matrix for each 'sequence.column' separately
      string_dist_matrices <- lapply(sequence_columns, function(x){
        
        # Extract sequences from 'nodes_list' of one 'sequence.column' and save in 'seqs' list
        seqs <- sapply(nodes_list, function(y) y[x])
        names(seqs) <- names(nodes_list)
        
        # Calculate distance matrix with 'stringdistmatrix()' function using the 'seqs' list as input 
        dist_matrix <- suppressWarnings(stringdist::stringdistmatrix(seqs, seqs, method = distance.metric, useNames = "names"))
        
        # Return the distance matrix
        return(dist_matrix)
      })
      
      # Rename distance matrices in list 'dist_matrices' according to the sequence columns in 'sequence_columns'
      names(string_dist_matrices) <- sequence_columns
    }
    
    
    # 3. Convert distance matrices to a lineage tree
    
    # If the provided 'network.algorithm' is a distance-based tree construction method, the function 'convert_dist_to_tree()' is used to build the lineage tree 
    if(network.algorithm %in% c("network.tree", "network.mst", "phylo.nj")){
      igraph_tree <- convert_dist_to_tree(clone = clone,
                                          dist.matrix = Reduce(`+`, string_dist_matrices), 
                                          network.algorithm = network.algorithm,
                                          resolve.ties = resolve.ties,
                                          node.sizes = node_sizes)
    }
    
    
    # 4. Create and return a list containing the phylogenetic tree, the sequences and features of the nodes in the tree, and the distance matrices used to build the tree
    
    # Return 'nodes_list', 'dist_matrices' and 'tree' as a list
    return(list(nodes = nodes_list, 
                distance.matrices = string_dist_matrices,
                lineage.tree = igraph_tree))
  }
  
  
  # Create list with sample IDs
  sample_list <- unique(VDJ$sample_id)
  
  # Create list with clonotype IDs 
  clone_list <- unique(paste(VDJ$sample_id, VDJ$clonotype_id, sep="_"))
  
  # Define partial function for to be executed for each clone
  partial_function <- function(clone){
    single_network <- infer_single_network(VDJ = VDJ,
                                           clone = clone,
                                           sequence.columns = sequence.columns,
                                           germline.columns = germline.columns,
                                           concatenate.sequences = concatenate.sequences,
                                           node.features = node.features,
                                           network.algorithm = network.algorithm,
                                           distance.metric = distance.metric,
                                           substitution.model = substitution.model,
                                           resolve.ties = resolve.ties)
    return(single_network)
  }
  
  # If 'parallel' is set to TRUE, the network inference is parallelized across the clones
  if(parallel){
    
    # Retrieve the operating system
    operating_system <- Sys.info()[['sysname']]
    
    # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
    if(operating_system %in% c('Linux', 'Darwin')){
      # Infer network for each clone 'clone_list' and store output in the list 'output_list'
      output_list <- parallel::mclapply(mc.cores = num.cores, X = clone_list, FUN = partial_function)
      # Name items in 'output_list' according to the clonotype IDs in 'clone_list'
      names(output_list) <- clone_list
    }
    
    # If the operating system is Windows, "parLapply" is used for parallelization
    if(operating_system == "Windows"){
      # Create cluster
      cluster <- parallel::makeCluster(num.cores)
      # Infer network for each clone 'clone_list' and store output in the list 'output_list'
      output_list <- parallel::parLapply(cluster, X = clone_list, fun = partial_function)
      # Name items in 'output_list' according to the clonotype IDs in 'clone_list'
      names(output_list) <- clone_list
      # Stop cluster
      parallel::stopCluster(cluster)
    }
    
  }
  
  # If 'parallel' is set to FALSE, the network inference is not parallelized
  if(!parallel){
    
    # Create a list for each sample and store in 'output_list'
    output_list <- lapply(clone_list, partial_function)
    # Name items in 'output_list' according to the clonotype IDs in 'clone_list'
    names(output_list) <- clone_list
    
  }
  
  # Reorganize list
  reorganized_output_list <- lapply(sample_list, function(sample){
    # Select all clonotypes from one sample
    matching_sublists <- grep(paste0("^",sample), names(output_list), value = TRUE)
    # Retrieve the networks from 'output_list' of one sample and store in 'sample_sublist'
    sample_sublist <- output_list[matching_sublists]
    # Rename the networks in 'sample_sublist' to their original clonotype ID
    names(sample_sublist) <- sub(paste0("^", sample, "_"), "", matching_sublists)
    # Return the 'sample_sublist'
    return(sample_sublist)
  })
  # Rename the sublists in 'reorganized_output_list' to their original sample ID
  names(reorganized_output_list) <- sample_list
  
  # Return 'output_list'
  return(reorganized_output_list)
}